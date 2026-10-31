#!/usr/bin/env nextflow
/*
 * RNASEQ - reproducible RNA-seq pipeline (paired-end; FASTQ pairs or BAM input), Docker version
 * (Windows / Mac personal computers). Same pipeline as the cluster version, run inside one Docker
 * image instead of LSF + Apptainer.
 * =============================================================================
 * [BAM -> BAM_TO_FASTQ] -> FastQC (raw) -> trim -> FastQC (trimmed) -> align (HISAT2) -> alignment summary
 * -> filter (MAPQ, remove_chroms) -> mark duplicates (or copy) -> gene model (GTF -> BED12) -> strandedness check
 * (RSeQC infer_experiment, warning only) -> featureCounts per sample -> merged counts, gene lengths, featureCounts
 * summary, samplesheet copy -> MultiQC.
 * Optional (run_deg: true): differential expression with DESeq2 on the merged counts table, or on an existing counts table
 * (deg_counts, which skips everything above): prepare -> QC plots -> per model DESeq2 fit -> results, volcano / MA plots,
 * heatmaps, gene plots.
 *
 * Run from your own project directory with run_rnaseq (run_rnaseq.ps1 or run_rnaseq.bat on Windows,
 * run_rnaseq.sh on Mac). It starts the pipeline's Docker image and runs nextflow inside it with the
 * -profile docker of nextflow.config:
 *   nextflow run /pipeline -profile docker -params-file config.yaml -resume
 * (docker/run_in_container.sh). Don't call nextflow yourself.
 *
 * Inputs (project root): config.yaml, samplesheet.csv
 */

nextflow.enable.dsl = 2

include { ALIGN } from './modules/local/align'
include { ALIGNMENT_SUMMARY } from './modules/local/alignment_summary'
include { BAM_FILTER } from './modules/local/bam_filter'
include { BAM_TO_FASTQ } from './modules/local/bam_to_fastq'
include { COUNTS_MERGE } from './modules/local/counts_merge'
include { DEG_DESEQ2 } from './modules/local/deg_deseq2'
include { DEG_GENE_PLOT } from './modules/local/deg_gene_plot'
include { DEG_HEATMAP } from './modules/local/deg_heatmap'
include { DEG_MA } from './modules/local/deg_ma'
include { DEG_PREPARE } from './modules/local/deg_prepare'
include { DEG_QC } from './modules/local/deg_qc'
include { DEG_RESULTS } from './modules/local/deg_results'
include { DEG_TRANSFORM } from './modules/local/deg_transform'
include { DEG_VOLCANO } from './modules/local/deg_volcano'
include { FASTQC as FASTQC_RAW } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { FEATURECOUNTS } from './modules/local/featurecounts'
include { GENEPRED_TO_BED12 } from './modules/local/gtf_to_bed12'
include { GTF_TO_GENEPRED } from './modules/local/gtf_to_bed12'
include { INFER_STRANDEDNESS } from './modules/local/infer_strandedness'
include { MARK_DUPLICATES } from './modules/local/mark_duplicates'
include { MULTIQC } from './modules/local/multiqc'
include { QC_STRANDEDNESS_TABLE } from './modules/local/qc_strandedness_table'
include { SAMPLESHEET_TSV } from './modules/local/samplesheet_tsv'
include { TRIM } from './modules/local/trim'

// -----------------------------------------------------------------------
// Parameter validation + CSV parsing, done once at run start.
// Nextflow's strict-syntax parser forbids bare top-level statements
// (if/for/etc.) outside a process, workflow or function, so all of this
// setup logic lives in functions called from the top of workflow{}.
// -----------------------------------------------------------------------

// Every step runs a script from this repo's bin/. A project replaces one by
// setting a config.yaml key named after the script (its file name without
// the extension) to the path of its own copy, e.g.
//     qc_frip: utils/qc_frip.sh
// Relative paths resolve against the project directory. Returns
// [script name: file] for every script in bin/, replacements applied.
def stepScripts() {
    def defaults = [:]
    def clashes = []
    file("${workflow.projectDir}/bin").listFiles().findAll { f -> f.isFile() && !f.name.startsWith('.') }.sort { f -> f.name }.each { f ->
        if (defaults.containsKey(f.baseName)) clashes << "${defaults[f.baseName].name} and ${f.name}"
        defaults[f.baseName] = f
    }
    if (clashes) fail("two files in bin/ share a base name (scripts are keyed by file name without extension): " + clashes.join('; '))

    // Keys from before scripts were replaced by name (none yet).
    def renamed = [:]
    def old = renamed.findAll { k, v -> params.containsKey(k) && params[k] }.collect { k, v -> "${k} -> ${v}" }
    if (old) fail("these config.yaml keys were renamed (a script is now replaced by a key named after it):\n  " + old.join('\n  '))

    // A script path under a key that names no pipeline script is most likely a typo.
    def stray = params.findAll { k, v -> v instanceof CharSequence && v.toString() =~ /(?i)\.(r|py|sh)$/ && !defaults.containsKey(k) }
    if (stray) {
        fail("config.yaml: no pipeline script is named ${stray.keySet().join(', ')}. To replace a step, name the key " +
             "after the script's file name without its extension. Scripts: ${defaults.keySet().sort().join(', ')}")
    }

    def scripts = [:]
    defaults.each { name, f ->
        // Only a string is a replacement path: params.trim is a boolean switch that shares its name with bin/trim.sh.
        def replacement = params.containsKey(name) && params[name] instanceof CharSequence ? params[name] : null
        if (replacement) {
            def repl = file(hostPathNote("script ${name}", replacement))
            if (!repl.exists()) fail("config.yaml: ${name}: ${replacement} not found.")
            log.info "Using this project's ${name}: ${replacement}"
            scripts[name] = repl
        } else {
            scripts[name] = f
        }
    }
    return scripts
}

def fail(String message) {
    log.error "ERROR: ${message}"
    System.exit(1)
}

// fail() with a title and one line per offending item; does nothing when items is empty.
def failList(String title, List items) {
    if (items) fail("${title}:\n  " + items.join('\n  '))
}

// Host paths. The launcher mounts the Windows drive X read-only at /mnt/x inside the Docker container, so a
// path written as D:\FASTQ\x.fastq.gz or D:/REF/hg.fa (config.yaml, CSV cells) is read as /mnt/d/FASTQ/x.fastq.gz.
// Any other string (cluster paths, relative paths, /Users/... on a Mac, UNC paths, blank) is returned unchanged.
// Wrap every user-supplied path with hostPath() where it becomes a file; error messages keep the path as written.
def hostPath(p) {
    if (p == null) return null
    def s = p.toString()
    def m = s =~ '^([A-Za-z]):[\\\\/](.*)$'
    if (!m.find()) return s
    return "/mnt/${m.group(1).toLowerCase()}/${m.group(2).replace('\\', '/')}".toString()
}

// The translations of this run: [where, as written, as used] rows. Kept in a JVM-wide list because strict
// syntax gives functions no access to script-level variables or `binding`; one Nextflow run per JVM.
def hostPathLog() {
    def rows = System.getProperties().get('rnaseq.host_paths')
    if (rows == null) {
        rows = []
        System.getProperties().put('rnaseq.host_paths', rows)
    }
    return rows
}

// hostPath() that also records a translated path (where = its param or CSV column) for reportHostPaths().
def hostPathNote(String where, p) {
    def used = hostPath(p)
    if (p != null && used != p.toString()) hostPathLog() << [where, p.toString(), used]
    return used
}

// Run start: when any Windows drive path was translated, logs the mapping and writes metadata/host_paths.tsv
// (key, as_written, as_used). Nothing is logged or written when no path was translated (the cluster).
def reportHostPaths() {
    def rows = hostPathLog().unique()
    if (!rows) return
    log.info "Windows drive paths are read from /mnt/<x> inside the container (host paths as written, and as used):"
    rows.each { r -> log.info "  ${r[0]}: ${r[1]} -> ${r[2]}" }
    try {
        def out = file("${params.metadata_dir}/host_paths.tsv")
        out.parent.mkdirs()
        out.text = (['key', 'as_written', 'as_used'].join('\t') + '\n') + rows.collect { r -> r.join('\t') }.join('\n') + '\n'
    } catch (Exception e) {
        log.warn "could not write ${params.metadata_dir}/host_paths.tsv: ${e.message}"
    }
}

// Container image of every process, as its module declares it
// (`container '<image>'`). Returns [process name: image].
def processContainers() {
    def images = [:]
    file("${projectDir}/modules/local").listFiles().findAll { f -> f.name.endsWith('.nf') }.sort { f -> f.name }.each { f ->
        def current = null
        f.readLines().each { line ->
            def p = line =~ '^\\s*process\\s+(\\w+)\\s*\\{'
            if (p.find()) current = p.group(1)
            def c = line =~ '^\\s*container\\s+[\'"]([^\'"$]+)[\'"]'
            if (current && c.find()) images[current] = c.group(1)
        }
    }
    return images
}

// true / 1 / yes (any case) or Boolean true; null and anything else is false.
def isTrue(v) {
    if (v == null) return false
    if (v instanceof Boolean) return v
    return v.toString().trim().toLowerCase() in ['true', '1', 'yes']
}

def isNum(v) {
    return v != null && v.toString().trim() ==~ /[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?/
}

// A positive integer (no sign, no decimal point; leading zeros allowed).
def isPosInt(v) {
    return textOf(v) ==~ /0*[1-9][0-9]*/
}

// Trimmed, non-empty items of a comma-separated value (null-safe).
def csvList(v) {
    if (v == null) return []
    return v.toString().split(',').collect { s -> s.trim() }.findAll { s -> s }
}

// Trimmed text of a config or CSV value ('' for null). Unlike ?:, keeps a numeric 0.
def textOf(v) {
    return v == null ? '' : v.toString().trim()
}

// Reads a small CSV: UTF-8 BOM and whitespace stripped from keys and values, fully blank rows
// dropped, every required column checked (values may be empty). Fails (exit 1) on a missing file
// or missing columns. Returns a list of row maps.
def readTable(String label, path, List required) {
    def f = file(hostPathNote(label, path))
    if (!f.exists()) fail("${label} not found: ${path}")
    def headerLine = f.readLines().find { l -> l.replace('﻿', '').trim() }
    def header = headerLine ? headerLine.split(',', -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() } : []
    def absent = required.findAll { c -> !(c in header) }
    if (absent) fail("${label} ${path} is missing column(s): ${absent.join(', ')}. Required columns: ${required.join(',')}")
    return f.splitCsv(header: true, quote: '"')
        .collect { row -> row.collectEntries { k, v -> [(k.toString().replace('﻿', '').trim()), (v ?: '').toString().trim()] } }
        .findAll { row -> row.values().any { v -> v } }
}

// HISAT2 index files of a prefix (<prefix>.N.ht2 / <prefix>.N.ht2l and the .rev. ones), sorted; empty when the index is missing.
def hisat2IndexFiles(String prefix) {
    def base = file(hostPath(prefix)).name
    return files("${hostPath(prefix)}*.ht2*").findAll { f -> f.name.startsWith("${base}.".toString()) }.sort { f -> f.name }
}

// config.yaml keys that are neither a declared parameter nor a script key (a bin/ script name, or <script>_args).
// The declared parameters are read from the params {} block of nextflow.config at run time (no second list to maintain).
def unknownKeys() {
    def declared = []
    def inParams = false
    file("${workflow.projectDir}/nextflow.config").readLines().each { line ->
        if (line ==~ /params\s*\{.*/) {
            inParams = true
        } else if (inParams && line ==~ /\}.*/) {
            inParams = false
        } else if (inParams) {
            def m = line =~ /^\s+(\w+)\s*=/
            if (m.find()) declared << m.group(1)
        }
    }
    def scripts = file("${workflow.projectDir}/bin").listFiles().findAll { f -> f.isFile() && !f.name.startsWith('.') }.collect { f -> f.baseName }
    def known = (declared + scripts + scripts.collect { s -> "${s}_args".toString() }) as Set
    // Keys with a '-' are Nextflow's own command-line spellings, not config.yaml keys.
    return params.keySet().collect { k -> k.toString() }.findAll { k -> !known.contains(k) && !k.contains('-') }.sort()
}

// An optional file input: the file, or an empty list when the param is unset (the process then stages nothing).
def optFile(v) {
    return textOf(v) ? file(hostPath(textOf(v))) : []
}

// true when v is empty or "a,b" with two numbers a < b (volcano axis limits).
def limitsOk(v) {
    def s = textOf(v)
    if (!s) return true
    def parts = s.split(',', -1).collect { p -> p.trim() }
    return parts.size() == 2 && parts.every { p -> isNum(p) } && parts[0].toDouble() < parts[1].toDouble()
}

// Header names of a small CSV / TSV (a .tsv or .txt file is tab separated), BOM and quotes removed; [] for an empty file.
def tableHeader(String path) {
    def f = file(hostPath(path))
    def line = f.readLines().find { l -> l.replace('﻿', '').trim() }
    def sep = f.name.toLowerCase().endsWith('.tsv') || f.name.toLowerCase().endsWith('.txt') ? '\t' : ','
    return line ? line.split(sep, -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() } : []
}

// Startup validation of a run_deg run: every deg_* param, the comparisons table and the heatmaps table. Everything is checked
// before any job starts. sheet / sheetHeader: the samplesheet path and header names when it is the sample metadata (deg_metadata
// unset), else null. Returns [metadata: path of the sample metadata, models: [[name, design, comparisons: [names]]] in table order,
// comparisons: the table rows, heatmaps: [row maps with every default filled in], geneX: gene plot column, geneList: [genes]].
def degConfig(sheet, sheetHeader) {
    def errs = []
    def idRe = '[A-Za-z0-9._+-]+'
    def boolOk = ['true', 'false']
    def hclustMethods = ['ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid']

    // ---- Switches and numbers ---------------------------------------------------
    ['deg_lfc_threshold_in_results', 'deg_transform_blind', 'deg_legacy_summary_counts', 'deg_run_qc', 'deg_run_barplot', 'deg_run_mds',
     'deg_run_correlation', 'deg_run_volcano', 'deg_run_volcano_pdfs', 'deg_run_volcano_combined', 'deg_run_heatmap', 'deg_run_ma',
     'deg_run_gene_plots'].each { k ->
        if (!(textOf(params[k]).toLowerCase() in boolOk)) errs << "${k}: '${textOf(params[k])}' must be true or false"
    }
    def alpha = isNum(params.deg_alpha) ? textOf(params.deg_alpha).toDouble() : -1d
    if (!(alpha > 0 && alpha < 1)) errs << "deg_alpha: '${textOf(params.deg_alpha)}' must be a number in (0, 1)"
    def lfc = isNum(params.deg_lfc) ? textOf(params.deg_lfc).toDouble() : -1d
    if (!(lfc >= 0)) errs << "deg_lfc: '${textOf(params.deg_lfc)}' must be a number >= 0"
    def sval = isNum(params.deg_sval_threshold) ? textOf(params.deg_sval_threshold).toDouble() : -1d
    if (!(sval > 0 && sval <= 1)) errs << "deg_sval_threshold: '${textOf(params.deg_sval_threshold)}' must be a number in (0, 1]"
    if (!isPosInt(params.deg_min_count)) errs << "deg_min_count: '${textOf(params.deg_min_count)}' must be a positive integer"
    if (!isPosInt(params.deg_min_samples)) errs << "deg_min_samples: '${textOf(params.deg_min_samples)}' must be a positive integer"
    [deg_shrink_type: ['ashr', 'apeglm', 'none'], deg_transform_method: ['rlog', 'vst'], deg_mds_input: ['counts', 'logcpm'],
     deg_cor_method: ['pearson', 'spearman', 'kendall'], deg_volcano_combined_layout: ['stacked', 'side'],
     deg_fit_type: ['parametric', 'local', 'mean']].each { k, ok ->
        if (!(textOf(params[k]) in ok)) errs << "${k}: '${textOf(params[k])}' is not valid (one of: ${ok.join(', ')})"
    }
    ['deg_volcano_xlim', 'deg_volcano_ylim', 'deg_volcano_common_xlim', 'deg_volcano_common_ylim'].each { k ->
        if (!limitsOk(params[k])) errs << "${k}: '${textOf(params[k])}' must be empty or two numbers 'min,max' with min < max"
    }
    ['deg_color_down', 'deg_color_up', 'deg_color_ns', 'deg_color_mid', 'deg_heatmap_color_up'].each { k ->
        if (!textOf(params[k])) errs << "${k} is empty (a colour name or a hex code such as '#8c6bb1')"
    }
    if (!textOf(params.deg_sample_col)) errs << "deg_sample_col is empty (the metadata column holding the sample names)"
    csvList(params.deg_reference_levels).each { item ->
        if (!(item ==~ /[^=]+=[^=]+/)) errs << "deg_reference_levels: '${item}' must be factor=level (comma separated list)"
    }

    // ---- Optional files -----------------------------------------------------------
    def counts = textOf(params.deg_counts)
    if (counts && !file(hostPathNote('deg_counts', counts)).exists()) errs << "deg_counts not found: ${counts}"
    ['deg_colors', 'deg_label_rules'].each { k ->
        if (textOf(params[k]) && !file(hostPathNote(k, textOf(params[k]))).exists()) errs << "${k} not found: ${textOf(params[k])}"
    }
    def metadata = sheet
    def metaOwn = textOf(params.deg_metadata)
    if (metaOwn) {
        metadata = metaOwn
        if (!file(hostPathNote('deg_metadata', metaOwn)).exists()) errs << "deg_metadata not found: ${metaOwn}"
        else {
            def cols = tableHeader(metaOwn)
            if (!(textOf(params.deg_sample_col) in cols))
                errs << "deg_metadata ${metaOwn} has no column '${textOf(params.deg_sample_col)}' (deg_sample_col); its columns: ${cols.join(', ')}"
        }
    } else if (!(textOf(params.deg_sample_col) in sheetHeader)) {
        errs << "samplesheet ${sheet} has no column '${textOf(params.deg_sample_col)}' (deg_sample_col); its columns: ${sheetHeader.join(', ')}"
    }

    // ---- Comparisons table --------------------------------------------------------
    if (!textOf(params.deg_comparisons)) fail("missing required param: deg_comparisons (run_deg is true: set the comparisons CSV in config.yaml).")
    def comps = readTable('deg_comparisons', textOf(params.deg_comparisons), ['model', 'design', 'comparison'])
    if (!comps) fail("deg_comparisons has no data rows: ${textOf(params.deg_comparisons)}")
    def models = []
    def seenKey = [] as Set
    comps.eachWithIndex { r, i ->
        def at = "deg_comparisons row ${i + 2} (${textOf(r.model)}/${textOf(r.comparison)})"   // CSV line number: the header is line 1
        if (!(textOf(r.model) ==~ idRe)) errs << "${at}: model '${textOf(r.model)}' must match ^${idRe}\$ (letters, digits, . _ + -; it names an output folder)"
        else if (textOf(r.model) in ['prepared', 'qc', 'gene_plots']) errs << "${at}: model '${textOf(r.model)}' is reserved (an output folder of its own); use another name"
        if (!(textOf(r.comparison) ==~ idRe)) errs << "${at}: comparison '${textOf(r.comparison)}' must match ^${idRe}\$ (letters, digits, . _ + -; it names an output folder and files)"
        else if (textOf(r.comparison) == 'log') errs << "${at}: comparison 'log' is reserved (the folder of the command logs)"
        if (!textOf(r.design).startsWith('~')) errs << "${at}: design '${textOf(r.design)}' must be a formula starting with ~ (e.g. ~treatment)"
        def hasChar = textOf(r.factor) || textOf(r.numerator) || textOf(r.denominator)
        def hasCoef = textOf(r.coef)
        def hasAdd = textOf(r.add) || textOf(r.subtract)
        def kinds = [hasChar, hasCoef, hasAdd].count { k -> k }
        if (kinds != 1) {
            errs << "${at}: needs exactly one contrast kind, found ${kinds > 1 ? 'more than one contrast kind' : 'none'} " +
                    "(factor+numerator+denominator, or coef, or add[+subtract])"
        } else if (hasChar && !(textOf(r.factor) && textOf(r.numerator) && textOf(r.denominator))) {
            errs << "${at}: a factor contrast needs all of factor, numerator and denominator"
        } else if (hasAdd && !textOf(r.add)) {
            errs << "${at}: subtract is set but add is empty (a list contrast needs add)"
        }
        if (textOf(params.deg_shrink_type) == 'apeglm' && !hasCoef)
            errs << "${at}: deg_shrink_type 'apeglm' needs a coef contrast (this row has none); use ashr or none, or a coef"
        if (textOf(r.combined_order) && !isNum(r.combined_order)) errs << "${at}: combined_order '${textOf(r.combined_order)}' must be a number"
        ['xlim', 'ylim'].each { k ->
            if (!limitsOk(r[k])) errs << "${at}: ${k} '${textOf(r[k])}' must be empty or two numbers 'min,max' with min < max"
        }
        if (!seenKey.add("${textOf(r.model)}\t${textOf(r.comparison)}".toString())) errs << "${at}: comparison '${textOf(r.comparison)}' appears twice in model '${textOf(r.model)}'"
    }
    comps.collect { r -> textOf(r.model) }.unique().each { m ->
        def designs = comps.findAll { r -> textOf(r.model) == m }.collect { r -> textOf(r.design) }.unique()
        if (designs.size() > 1) errs << "deg_comparisons: model '${m}' has more than one design (${designs.join(' | ')}): one model is one DESeq2 fit, give each design its own model name"
        models << [name: m, design: designs[0], comparisons: comps.findAll { r -> textOf(r.model) == m }.collect { r -> textOf(r.comparison) }]
    }

    // ---- Heatmaps table -------------------------------------------------------------
    def heatmaps = []
    if (textOf(params.deg_heatmaps)) {
        def hrows = readTable('deg_heatmaps', textOf(params.deg_heatmaps), ['name', 'model', 'comparison'])
        def names = hrows.collect { r -> textOf(r.name) }
        def dupNames = names.findAll { n -> n && names.count(n) > 1 }.unique()
        if (dupNames) errs << "deg_heatmaps: duplicate name(s): ${dupNames.join(', ')}"
        hrows.eachWithIndex { r, i ->
            def at = "deg_heatmaps row ${i + 2} (${textOf(r.name)})"
            if (!(textOf(r.name) ==~ idRe)) errs << "${at}: name '${textOf(r.name)}' must match ^${idRe}\$ (it names the output file)"
            def crow = comps.find { c -> textOf(c.model) == textOf(r.model) && textOf(c.comparison) == textOf(r.comparison) }
            if (!crow) errs << "${at}: no comparison '${textOf(r.comparison)}' in model '${textOf(r.model)}' of deg_comparisons"
            if (textOf(r.gene_set) && !file(hostPathNote('deg_heatmaps gene_set', textOf(r.gene_set))).exists()) errs << "${at}: gene_set not found: ${textOf(r.gene_set)}"
            ['cutree_rows', 'palette_n'].each { k ->
                if (textOf(r[k]) && !isPosInt(r[k])) errs << "${at}: ${k} '${textOf(r[k])}' must be a positive integer"
            }
            ['fontsize', 'fontsize_row', 'width', 'height'].each { k ->
                if (textOf(r[k]) && !(isNum(r[k]) && textOf(r[k]).toDouble() > 0)) errs << "${at}: ${k} '${textOf(r[k])}' must be a positive number"
            }
            if (textOf(r.clust_method) && !(textOf(r.clust_method) in hclustMethods))
                errs << "${at}: clust_method '${textOf(r.clust_method)}' is not valid (one of: ${hclustMethods.join(', ')})"
            if (textOf(r.col_swap) && !(textOf(r.col_swap) in ['none', 'scaled', 'unscaled']))
                errs << "${at}: col_swap '${textOf(r.col_swap)}' is not valid (one of: none, scaled, unscaled)"
            def label = crow && textOf(crow.label) ? textOf(crow.label) : textOf(r.comparison).replace('_', ' ')
            heatmaps << [name: textOf(r.name), model: textOf(r.model), comparison: textOf(r.comparison), gene_set: textOf(r.gene_set),
                         clust_method: textOf(r.clust_method) ?: 'ward.D2', cutree_rows: textOf(r.cutree_rows), col_swap: textOf(r.col_swap) ?: 'none',
                         palette_n: textOf(r.palette_n) ?: '1000', fontsize: textOf(r.fontsize) ?: '9', fontsize_row: textOf(r.fontsize_row) ?: '6.5',
                         width: textOf(r.width) ?: '6.5', height: textOf(r.height) ?: '4.5',
                         title: (textOf(r.title) ?: "DEGs: ${label} (N = {n})").toString()]
        }
        // Two different gene set files with one file name would be confusing in the logs and the work directories.
        def gsPaths = heatmaps.collect { h -> h.gene_set }.findAll { p -> p }.unique()
        gsPaths.collect { p -> file(hostPath(p)).name }.unique().each { n ->
            def same = gsPaths.findAll { p -> file(hostPath(p)).name == n }
            if (same.size() > 1) errs << "deg_heatmaps: gene_set files must have unique file names; '${n}' is used by: ${same.join(', ')}"
        }
    }

    failList("invalid differential expression (deg_*) setting(s)", errs)

    // Gene plots: the column on the x axis and for the colours.
    def designVars = models[0].design.replaceFirst(/^\s*~/, '').split(/[^A-Za-z0-9._]+/).findAll { v -> v }
    def geneX = textOf(params.deg_gene_plot_x) ?: textOf(params.deg_color_by) ?: (designVars ? designVars[0] : '')
    return [metadata: metadata, models: models, comparisons: comps, heatmaps: heatmaps, geneX: geneX, geneList: csvList(params.deg_gene_plot_genes)]
}

// DEG lines of the startup summary box (empty when run_deg is false).
def degSummary(deg, boolean countsOnly) {
    if (!deg) return ''
    def src = countsOnly ? "counts table ${params.deg_counts}" : 'merged counts of this run'
    def nComp = deg.comparisons.size()
    def shrink = textOf(params.deg_shrink_type)
    def parts = []
    if (isTrue(params.deg_run_qc)) parts << 'qc'
    if (isTrue(params.deg_run_volcano)) parts << 'volcano'
    if (isTrue(params.deg_run_ma)) parts << 'ma'
    if (isTrue(params.deg_run_heatmap) && deg.heatmaps) parts << "heatmaps(${deg.heatmaps.size()})"
    if (isTrue(params.deg_run_gene_plots) && deg.geneList) parts << "gene plots(${deg.geneList.size()})"
    return "\n" +
        "      DEG            : ${nComp} comparison(s) in ${deg.models.size()} model(s) (${deg.models.collect { m -> m.name }.join(', ')}); input ${src}\n" +
        "      DEG settings   : alpha ${textOf(params.deg_alpha)}, |log2FC| ${textOf(params.deg_lfc)}, shrinkage ${shrink}, ${textOf(params.deg_transform_method)} for heatmaps\n" +
        "      DEG outputs    : ${parts.join(', ') ?: 'results tables only'}"
}

def prepareRun() {
    if (params.help) {
        log.info """
        RNASEQ pipeline (RNA-seq, paired-end; FASTQ or BAM input)
        ---------------------------------------------------------
        Usage:
          nextflow run ${workflow.projectDir} -profile docker -params-file config.yaml

        Required params (set in config.yaml):
          samplesheet       CSV: sample_id, fastq_1 + fastq_2 or bam, exclude, plus free metadata columns
                            (sample_id required; paired-end only; one input kind per row)
          hisat2_index      HISAT2 index PREFIX of the genome (files <prefix>.N.ht2 or <prefix>.N.ht2l)
          gtf               gene annotation GTF

        Optional params:
          strandedness      unstranded (default) | forward | reverse (featureCounts library type; check reported in QC)
          run_fastqc, trim, mark_duplicates, remove_duplicates, run_multiqc: true | false
          mapq, remove_chroms, hisat2_args, featurecounts_args, strand_call_threshold, ...
          project_name, outdir (./results)

        Differential expression (optional, DESeq2; off by default):
          run_deg           true | false. true adds DEG_* after the counts table
          deg_comparisons   CSV, one row per comparison (required with run_deg): model, design, comparison + contrast columns
          deg_counts        an existing counts table: with run_deg it skips ALL upstream steps (no reads, QC or MultiQC)
          deg_metadata      sample metadata CSV/TSV (default: the samplesheet); deg_heatmaps, deg_colors, deg_label_rules: optional CSVs
          deg_alpha, deg_lfc, deg_shrink_type, deg_transform_method, deg_run_* switches, ...: README.md ("Differential expression")

        The other keys and their defaults are in nextflow.config. Replace one step's script with
        <script name>: path in config.yaml (see stepScripts()); extra arguments with <script name>_args.
        """
        System.exit(0)
    }

    // Keys that nothing reads (a typo would otherwise be silently ignored).
    def unknown = unknownKeys()
    if (unknown) fail("unknown config key(s): ${unknown.join(', ')} (not a pipeline parameter; not a script key or <script>_args; check the spelling against the params block of nextflow.config)")

    def errs = []
    def idRe = '[A-Za-z0-9._-]+'
    def boolOk = ['true', 'false']

    // ---- Top-level params ---------------------------------------------------
    ['run_fastqc', 'trim', 'run_multiqc', 'mark_duplicates', 'remove_duplicates', 'run_deg'].each { k ->
        if (!(textOf(params[k]).toLowerCase() in boolOk)) errs << "${k}: '${textOf(params[k])}' must be true or false"
    }
    if (!(textOf(params.strandedness) in ['unstranded', 'forward', 'reverse']))
        errs << "strandedness: '${textOf(params.strandedness)}' is not valid (one of: unstranded, forward, reverse)"
    if (!(textOf(params.mapq) ==~ /[0-9]+/)) errs << "mapq: '${textOf(params.mapq)}' must be a non-negative integer"
    if (!isPosInt(params.trim_minlen)) errs << "trim_minlen: '${textOf(params.trim_minlen)}' must be a positive integer"
    if (!isPosInt(params.optical_duplicate_pixel_distance))
        errs << "optical_duplicate_pixel_distance: '${textOf(params.optical_duplicate_pixel_distance)}' must be a positive integer"
    if (!isPosInt(params.infer_experiment_reads)) errs << "infer_experiment_reads: '${textOf(params.infer_experiment_reads)}' must be a positive integer"
    def alignPct = isNum(params.alignment_min_pct) ? textOf(params.alignment_min_pct).toDouble() : -1d
    if (!(alignPct >= 0 && alignPct <= 100)) errs << "alignment_min_pct: '${textOf(params.alignment_min_pct)}' must be a number between 0 and 100"
    def strandThr = isNum(params.strand_call_threshold) ? textOf(params.strand_call_threshold).toDouble() : -1d
    if (!(strandThr > 0.5 && strandThr <= 1)) errs << "strand_call_threshold: '${textOf(params.strand_call_threshold)}' must be a number in (0.5, 1]"
    if (!textOf(params.featurecounts_feature_type)) errs << "featurecounts_feature_type is empty (e.g. exon)"
    if (!textOf(params.featurecounts_attribute)) errs << "featurecounts_attribute is empty (e.g. gene_id)"
    if (isTrue(params.remove_duplicates) && !isTrue(params.mark_duplicates))
        errs << "remove_duplicates requires mark_duplicates (remove_duplicates is true but mark_duplicates is false)"
    failList("invalid config.yaml value(s)", errs)

    // ---- Differential expression mode -------------------------------------------------
    // deg_counts (with run_deg) runs the DEG steps on an existing counts table: no reads, references, QC or MultiQC are needed.
    def degOn = isTrue(params.run_deg)
    def countsOnly = degOn && textOf(params.deg_counts) ? true : false
    if (textOf(params.deg_counts) && !degOn) fail("deg_counts is set but run_deg is false (deg_counts only feeds the differential expression: set run_deg: true).")
    if (countsOnly && textOf(params.deg_metadata)) {
        // The sample metadata comes from deg_metadata: the samplesheet is not read at all.
        def degOnly = degConfig(null, null)
        log.info """
    RNASEQ  ~  ${params.project_name}  (differential expression only)
    -----------------------------------------
      counts table   : ${params.deg_counts}
      metadata       : ${params.deg_metadata}
      outdir         : ${params.outdir}${degSummary(degOnly, true)}
    """
        return [scripts: stepScripts(), images: processContainers(), deg: degOnly, countsOnly: true]
    }

    if (!textOf(params.samplesheet)) fail("missing required param: samplesheet (set in your config.yaml).")

    // ---- Samplesheet ----------------------------------------------------------
    def sheet = params.samplesheet.toString()
    def sheetFile = file(hostPathNote('samplesheet', sheet))
    if (!sheetFile.exists()) fail("samplesheet not found: ${sheet}")
    // Checked on the raw header first: a duplicated column name would break the CSV parser.
    def headerLine = sheetFile.readLines().find { l -> l.replace('﻿', '').trim() }
    def header = headerLine ? headerLine.split(',', -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() } : []
    def dupCols = header.findAll { c -> c && header.count(c) > 1 }.unique()
    if (dupCols) fail("samplesheet ${sheet}: duplicate column header(s): ${dupCols.join(', ')}")
    def allRows = readTable('samplesheet', sheet, ['sample_id'])
    if (!allRows) fail("samplesheet has no data rows: ${sheet}")
    // Optional input and flag columns: blank when absent. Other columns are free metadata.
    allRows.each { row -> ['fastq_1', 'fastq_2', 'bam', 'exclude'].each { c -> row[c] = row[c] ?: '' } }

    // sample_id: present, safe in file names, unique.
    allRows.eachWithIndex { row, i ->
        if (!row.sample_id) errs << "row ${i + 1} has an empty sample_id"
        else if (!(row.sample_id ==~ idRe)) errs << "sample_id '${row.sample_id}' must match ^${idRe}\$ (letters, digits, _ . -; used verbatim in file names)"
    }
    def allIds = allRows.collect { row -> row.sample_id }
    def dupIds = allIds.findAll { s -> s && allIds.count(s) > 1 }.unique()
    if (dupIds) errs << "duplicate sample_id values: ${dupIds.join(', ')}"

    // exclude flag.
    allRows.each { row ->
        def e = row.exclude.toLowerCase()
        if (e in ['', 'false', '0', 'no']) row.excluded = false
        else if (e in ['true', '1', 'yes']) row.excluded = true
        else {
            errs << "${row.sample_id}: exclude '${row.exclude}' must be empty/false/0/no or true/1/yes"
            row.excluded = false
        }
    }

    // Only the rows that are not excluded are processed.
    def rows = allRows.findAll { row -> !row.excluded }
    def excludedIds = allRows.findAll { row -> row.excluded }.collect { row -> row.sample_id }
    if (!rows) fail("every row is excluded: ${sheet}")
    def ids = rows.collect { row -> row.sample_id }

    // Input per row: one FASTQ pair (paired-end only) or one BAM (not checked when only the counts table is used: deg_counts).
    (countsOnly ? [] : rows).each { row ->
        def id = row.sample_id
        def hasFq = (row.fastq_1 || row.fastq_2) ? true : false
        def hasBam = row.bam ? true : false
        if (hasFq && hasBam) errs << "${id}: both fastq_1/fastq_2 and bam are set (give a FASTQ pair or a BAM, not both)"
        else if (!hasFq && !hasBam) errs << "${id}: neither fastq_1/fastq_2 nor bam is set"
        else if (hasBam) {
            if (!file(hostPathNote("samplesheet ${id} bam", row.bam)).exists()) errs << "${id}: bam not found: ${row.bam}"
        } else if (!row.fastq_1 || !row.fastq_2) {
            errs << "${id}: single-end input is not supported (paired-end only): fastq_1 and fastq_2 are both required"
        } else {
            ['fastq_1', 'fastq_2'].each { k -> if (!file(hostPathNote("samplesheet ${id} ${k}", row[k])).exists()) errs << "${id}: ${k} not found: ${row[k]}" }
            if (row.fastq_1 == row.fastq_2) errs << "${id}: fastq_1 and fastq_2 are the same file (${row.fastq_1})"
        }
    }
    def fqPaths = (countsOnly ? [] : rows).findAll { row -> !row.bam }.collectMany { row -> [row.fastq_1, row.fastq_2] }.findAll { p -> p }
    def fqNames = fqPaths.collect { p -> file(hostPath(p)).name }
    fqNames.findAll { n -> fqNames.count(n) > 1 }.unique().each { n ->
        errs << "two FASTQ files share the name ${n}: " + fqPaths.findAll { p -> file(hostPath(p)).name == n }.join(', ')
    }
    failList("${sheet}: invalid sample / input values", errs)

    // Differential expression: the deg_* params and the comparisons / heatmaps tables (the samplesheet is the sample metadata unless deg_metadata is set).
    def deg = degOn ? degConfig(sheet, header) : null

    def fastqRows = rows.findAll { row -> !row.bam }
    def bamRows = rows.findAll { row -> row.bam }

    // ---- Reference files ------------------------------------------------------
    def refs = [:]
    if (!countsOnly) {
        def refErrs = []
        def prefix = textOf(params.hisat2_index)
        def indexFiles = []
        def indexBase = null
        if (!prefix) {
            refErrs << "hisat2_index: required (HISAT2 index prefix, e.g. /path/to/t2t/chm13-2)"
        } else {
            indexBase = file(hostPathNote('hisat2_index', prefix)).name
            indexFiles = hisat2IndexFiles(prefix)
            if (!indexFiles.any { f -> f.name == "${indexBase}.1.ht2" || f.name == "${indexBase}.1.ht2l" })
                refErrs << "hisat2_index: no HISAT2 index files (${prefix}.1.ht2 or .1.ht2l) found for prefix ${prefix}"
        }
        def gtf = textOf(params.gtf)
        if (!gtf) refErrs << "gtf: required (gene annotation GTF)"
        else if (!file(hostPathNote('gtf', gtf)).exists()) refErrs << "gtf not found: ${gtf}"
        failList("reference file problems", refErrs)
        refs = [index_files: indexFiles, index_base: indexBase, gtf: file(hostPath(gtf))]
    }

    if (countsOnly) {
        // run_deg with deg_counts and no deg_metadata: the samplesheet only supplies the sample metadata.
        log.info """
    RNASEQ  ~  ${params.project_name}  (differential expression only)
    -----------------------------------------
      counts table   : ${params.deg_counts}
      metadata       : samplesheet ${sheet}: ${ids.size()} active row(s), ${excludedIds.size()} excluded
      outdir         : ${params.outdir}${degSummary(deg, true)}
    """
        return [scripts: stepScripts(), images: processContainers(), deg: deg, countsOnly: true]
    }

    // ---- Summary --------------------------------------------------------------
    def steps = []
    ['run_fastqc', 'trim', 'mark_duplicates', 'remove_duplicates', 'run_multiqc'].each { k ->
        if (isTrue(params[k])) steps << k
    }
    def dedupMode = 'mark only'
    if (!isTrue(params.mark_duplicates)) dedupMode = 'off (the filtered BAM is the final BAM)'
    else if (isTrue(params.remove_duplicates)) dedupMode = 'mark and remove'
    log.info """
    RNASEQ  ~  ${params.project_name}
    -----------------------------------------
      samples        : ${ids.size()} active (FASTQ ${fastqRows.size()}, BAM ${bamRows.size()}), ${excludedIds.size()} excluded${excludedIds ? ' (' + excludedIds.join(', ') + ')' : ''}
      strandedness   : ${textOf(params.strandedness)} (strandedness check threshold ${textOf(params.strand_call_threshold)})
      MAPQ           : ${textOf(params.mapq)}
      duplicates     : ${dedupMode}
      steps on       : ${steps.join(', ')}
      outdir         : ${params.outdir}${degSummary(deg, countsOnly)}
    """

    return [scripts: stepScripts(), images: processContainers(), refs: refs, fastqRows: fastqRows, bamRows: bamRows,
            excludedIds: excludedIds, ids: ids, deg: deg, countsOnly: countsOnly]
}

workflow {

    def cfg = prepareRun()
    reportHostPaths()

    // Each step's script: bin/<name>, or this project's replacement (see stepScripts()).
    def sc = cfg.scripts
    // The counts table of the DEG steps: this run's merged counts, or params.deg_counts.
    def counts_table = null

    // Reads, QC, alignment and counting. Skipped when only the counts table is used (run_deg with deg_counts).
    if (!cfg.countsOnly) {
        def samplesheet = Channel.value(file(hostPath(params.samplesheet)))
        def gtf = Channel.value(cfg.refs.gtf)
        // featureCounts -s code for each library type.
        def strandCode = ['unstranded': 0, 'forward': 1, 'reverse': 2]

        // ---- 0. Reads: FASTQ pairs as given; BAMs converted to FASTQ pairs --------------
        def reads = Channel.fromList(cfg.fastqRows.collect { row -> tuple(row.sample_id, file(hostPath(row.fastq_1)), file(hostPath(row.fastq_2))) })
        if (cfg.bamRows) {
            BAM_TO_FASTQ(Channel.fromList(cfg.bamRows.collect { row -> tuple(row.sample_id, file(hostPath(row.bam))) }), sc.bam_to_fastq)
            reads = reads.mix(BAM_TO_FASTQ.out.reads)
        }

        // ---- 1. FastQC (raw), trimming, FastQC (trimmed) ---------------------------------
        def fastqc_zips = Channel.empty()
        def trim_logs = Channel.empty()
        if (isTrue(params.run_fastqc)) {
            FASTQC_RAW(reads.flatMap { id, r1, r2 -> [tuple(id, r1), tuple(id, r2)] }, 'raw', sc.fastqc)
            fastqc_zips = FASTQC_RAW.out.zip.map { id, zip -> zip }
        }
        if (isTrue(params.trim)) {
            TRIM(reads, params.trim_adapter_set, params.trim_illuminaclip, params.trim_steps, params.trim_minlen, sc.trim)
            reads = TRIM.out.reads
            trim_logs = TRIM.out.log
            if (isTrue(params.run_fastqc)) {
                FASTQC_TRIM(reads.flatMap { id, r1, r2 -> [tuple(id, r1), tuple(id, r2)] }, 'trimmed', sc.fastqc)
                fastqc_zips = fastqc_zips.mix(FASTQC_TRIM.out.zip.map { id, zip -> zip })
            }
        }

        // ---- 2. Align (HISAT2 | samtools), alignment summary, filter, duplicates ------------
        ALIGN(reads, Channel.value(cfg.refs.index_files), cfg.refs.index_base, textOf(params.hisat2_args), sc.align)
        ALIGNMENT_SUMMARY(ALIGN.out.summary.collect(), params.alignment_min_pct, sc.alignment_summary)
        BAM_FILTER(ALIGN.out.bam, params.mapq, textOf(params.remove_chroms), sc.bam_filter)
        // Always run: with mark_duplicates false the task only copies the filtered BAM to S.final.bam.
        MARK_DUPLICATES(BAM_FILTER.out.bam, isTrue(params.mark_duplicates).toString(), isTrue(params.remove_duplicates).toString(),
                        params.optical_duplicate_pixel_distance, sc.mark_duplicates)

        // ---- 3. Strandedness: gene model, per-sample check, one table (warning only) ----------
        GTF_TO_GENEPRED(gtf, sc.gtf_to_genepred)
        GENEPRED_TO_BED12(GTF_TO_GENEPRED.out.genepred, sc.genepred_to_bed12)
        INFER_STRANDEDNESS(MARK_DUPLICATES.out.bam, GENEPRED_TO_BED12.out.bed12, params.infer_experiment_reads, sc.infer_strandedness)
        QC_STRANDEDNESS_TABLE(INFER_STRANDEDNESS.out.txt.collect(), textOf(params.strandedness), params.strand_call_threshold, sc.qc_strandedness_table)

        // ---- 4. Counts: featureCounts per sample, merged table, samplesheet copy ------------
        FEATURECOUNTS(MARK_DUPLICATES.out.bam, gtf, strandCode[textOf(params.strandedness)], params.featurecounts_feature_type,
                      params.featurecounts_attribute, sc.featurecounts)
        def per_sample = FEATURECOUNTS.out.counts.mix(FEATURECOUNTS.out.table, FEATURECOUNTS.out.summary).collect()
        COUNTS_MERGE(per_sample, cfg.ids.join(','), params.project_name, sc.counts_merge)
        SAMPLESHEET_TSV(samplesheet, cfg.ids.join(','), sc.samplesheet_tsv)

        // ---- 5. MultiQC ------------------------------------------------------------------
        if (isTrue(params.run_multiqc)) {
            def multiqc_in = fastqc_zips.mix(trim_logs, ALIGN.out.summary, BAM_FILTER.out.raw_flagstat, BAM_FILTER.out.raw_idxstats,
                                             BAM_FILTER.out.filtered_flagstat, BAM_FILTER.out.filtered_idxstats,
                                             MARK_DUPLICATES.out.metrics, FEATURECOUNTS.out.summary, INFER_STRANDEDNESS.out.txt).collect()
            MULTIQC(multiqc_in, sc.multiqc)
        }

        counts_table = COUNTS_MERGE.out.counts
    }

    // ---- 6. Differential expression (run_deg) -----------------------------------------------
    // DEG_PREPARE once (counts, metadata and comparisons tidied), DEG_QC once, then per model of deg_comparisons: DEG_DESEQ2 ->
    // DEG_RESULTS -> DEG_VOLCANO / DEG_MA. DEG_TRANSFORM (rlog / vst) runs only for the models that heatmaps or gene plots need;
    // DEG_HEATMAP runs once per row of deg_heatmaps, DEG_GENE_PLOT once (on the first model's transformed counts).
    if (isTrue(params.run_deg)) {
        def deg = cfg.deg
        if (cfg.countsOnly) counts_table = Channel.value(file(hostPath(params.deg_counts)))
        def metadata_in = Channel.value(file(hostPath(deg.metadata)))
        def comparisons_in = Channel.value(file(hostPath(params.deg_comparisons)))
        // Optional files: the file, or an empty list (the module then stages nothing).
        def colors_in = Channel.value(optFile(params.deg_colors))
        def label_rules_in = Channel.value(optFile(params.deg_label_rules))

        DEG_PREPARE(counts_table, metadata_in, comparisons_in, sc.deg_common, textOf(params.deg_sample_col), textOf(params.deg_strip_suffix),
                    textOf(params.deg_reference_levels), sc.deg_prepare)
        // DEG_PREPARE's inputs are all value channels, so its outputs are value channels too: read by one task per model (and per heatmap).
        def prep_counts = DEG_PREPARE.out.counts
        def prep_metadata = DEG_PREPARE.out.metadata
        def prep_comparisons = DEG_PREPARE.out.comparisons

        if (isTrue(params.deg_run_qc)) {
            DEG_QC(prep_counts, prep_metadata, prep_comparisons, colors_in, sc.deg_common, textOf(params.deg_color_by),
                   isTrue(params.deg_run_barplot).toString(), isTrue(params.deg_run_mds).toString(), isTrue(params.deg_run_correlation).toString(),
                   textOf(params.deg_mds_input), textOf(params.deg_mds_label_by), textOf(params.deg_mds_legend_title),
                   textOf(params.deg_legacy_mds_legend_col), textOf(params.deg_cor_method), sc.deg_qc)
        }

        // One DESeq2 fit per model (rows of deg_comparisons with the same model share it), then its results tables.
        def models = Channel.fromList(deg.models.collect { m -> tuple(m.name, m.design) })
        DEG_DESEQ2(models, prep_counts, prep_metadata, sc.deg_common, textOf(params.deg_min_count), textOf(params.deg_min_samples),
                   textOf(params.deg_fit_type), sc.deg_deseq2)
        DEG_RESULTS(DEG_DESEQ2.out.deseq, prep_comparisons, sc.deg_common, textOf(params.deg_alpha), textOf(params.deg_lfc),
                    isTrue(params.deg_lfc_threshold_in_results).toString(), textOf(params.deg_shrink_type), textOf(params.deg_sval_threshold),
                    isTrue(params.deg_legacy_summary_counts).toString(), sc.deg_results)

        if (isTrue(params.deg_run_volcano)) {
            DEG_VOLCANO(DEG_RESULTS.out.results, prep_comparisons, sc.deg_common, textOf(params.deg_alpha), textOf(params.deg_lfc),
                        textOf(params.deg_sval_threshold), (textOf(params.deg_shrink_type) != 'none').toString(), textOf(params.deg_color_down),
                        textOf(params.deg_color_up), textOf(params.deg_color_ns), textOf(params.deg_label_genes), textOf(params.deg_volcano_xlim),
                        textOf(params.deg_volcano_ylim), textOf(params.deg_volcano_common_xlim), textOf(params.deg_volcano_common_ylim),
                        textOf(params.deg_volcano_combined_layout), isTrue(params.deg_run_volcano_pdfs).toString(),
                        isTrue(params.deg_run_volcano_combined).toString(), sc.deg_volcano)
        }
        if (isTrue(params.deg_run_ma)) {
            DEG_MA(DEG_RESULTS.out.results, prep_comparisons, sc.deg_common, textOf(params.deg_alpha), textOf(params.deg_lfc),
                   textOf(params.deg_color_down), textOf(params.deg_color_up), textOf(params.deg_color_ns), sc.deg_ma)
        }

        // Transformed counts (slow: rlog) only for the models that need them: the models of the heatmaps, and the first model for gene plots.
        def heatmaps_on = isTrue(params.deg_run_heatmap) && deg.heatmaps
        def gene_plots_on = isTrue(params.deg_run_gene_plots) && deg.geneList
        def heatmap_models = heatmaps_on ? deg.heatmaps.collect { h -> h.model }.unique() : []
        def transform_models = heatmap_models + (gene_plots_on ? [deg.models[0].name] : [])
        if (transform_models) {
            DEG_TRANSFORM(DEG_DESEQ2.out.dataset.filter { m, f -> m in transform_models }, sc.deg_common, textOf(params.deg_transform_method),
                          isTrue(params.deg_transform_blind).toString(), sc.deg_transform)

            // One heatmap per row of deg_heatmaps: the comparison's results table and the model's transformed counts, keyed by model.
            if (heatmaps_on) {
                def results_for_heatmaps = DEG_RESULTS.out.results.filter { m, f -> m in heatmap_models }
                def counts_for_heatmaps = DEG_TRANSFORM.out.counts.filter { m, f -> m in heatmap_models }
                def heatmap_inputs = results_for_heatmaps.join(counts_for_heatmaps, failOnMismatch: true, failOnDuplicate: true)
                def heatmap_rows = Channel.fromList(deg.heatmaps.collect { h -> tuple(h.model, h) })
                    .combine(heatmap_inputs, by: 0)
                    .map { model, h, results, rcounts ->
                        def one = [results].flatten().find { f -> f.name == "results_${h.comparison}.csv" }
                        tuple(h.findAll { k, v -> k != 'gene_set' }, one, rcounts, h.gene_set ? file(hostPath(h.gene_set)) : [])
                    }
                DEG_HEATMAP(heatmap_rows, label_rules_in, sc.deg_common, textOf(params.deg_color_down), textOf(params.deg_heatmap_color_up),
                            textOf(params.deg_color_mid), sc.deg_heatmap)
            }

            // Expression plots of single genes, from the first model's transformed counts.
            if (gene_plots_on) {
                def first_model = deg.models[0].name
                DEG_GENE_PLOT(DEG_TRANSFORM.out.counts.filter { m, f -> m == first_model }.map { m, f -> f }, prep_metadata, colors_in,
                              sc.deg_common, deg.geneList.join(','), deg.geneX, sc.deg_gene_plot)
            }
        }
    }
}
