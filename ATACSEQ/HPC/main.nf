#!/usr/bin/env nextflow
/*
 * ATACSEQ - reproducible ATAC-seq pipeline (paired-end)
 * ======================================================
 * FASTQ -> FastQC -> trim -> align (bowtie2, optionally in read-pair chunks merged per library) -> alignment summary -> filter (MAPQ, remove_chroms) -> dedup ->
 * blacklist -> final BAM -> [contamination] -> [Tn5 shift] -> QC (fragment sizes, chrM fraction) -> TMM scale
 * factors -> bigWigs -> TSS enrichment -> MACS3 per peak group (+ peak merges) -> FRiP, peak stats -> [differential accessibility,
 * peak annotation, motifs: run_differential, run_annotation, run_motifs] -> MultiQC.
 * A samplesheet row may carry a final BAM (column bam) instead of FASTQ; it joins the final BAMs of the FASTQ libraries (BAM_PREALIGNED).
 *
 * Run from your own project directory with:
 *   nextflow run /project/rephilllab/ATACSEQ -profile pmacs -params-file config.yaml
 *
 * Inputs (project root): config.yaml, samplesheet.csv, peak_merges.csv (optional), and for the downstream steps comparisons.csv,
 * peak_files.csv, diff_tables.csv, motif_targets.csv (each optional, named in config.yaml)
 */

nextflow.enable.dsl = 2

include { ALIGN } from './modules/local/align'
include { ALIGNMENT_SUMMARY } from './modules/local/alignment_summary'
include { ALIGN_MERGE } from './modules/local/align_merge'
include { BAM_FILTER } from './modules/local/bam_filter'
include { BAM_FINALIZE } from './modules/local/bam_finalize'
include { BAM_PREALIGNED } from './modules/local/bam_prealigned'
include { BIGWIG } from './modules/local/bigwig'
include { BLACKLIST_FILTER } from './modules/local/blacklist_filter'
include { CONTAMINATION } from './modules/local/contamination'
include { COUNT_PAIRS } from './modules/local/count_pairs'
include { DEDUP } from './modules/local/dedup'
include { DIFF_ACCESS } from './modules/local/diff_access'
include { DIFF_ANNOTATE } from './modules/local/diff_annotate'
include { FASTQC as FASTQC_RAW } from './modules/local/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/local/fastqc'
include { FRIP } from './modules/local/frip'
include { MACS3 } from './modules/local/macs3'
include { MOTIF_COMMON } from './modules/local/motif_common'
include { MOTIF_ENRICHMENT } from './modules/local/motif_enrichment'
include { MULTIQC } from './modules/local/multiqc'
include { NORM_TMM } from './modules/local/norm_tmm'
include { PEAK_MERGE } from './modules/local/peak_merge'
include { QC_CHROM_FRACTION } from './modules/local/qc_chrom_fraction'
include { QC_FRAGMENT_SIZES_PLOT } from './modules/local/qc_fragment_sizes_plot'
include { QC_FRAGMENT_SIZES_TABLE } from './modules/local/qc_fragment_sizes_table'
include { QC_FRIP_TABLE } from './modules/local/qc_frip_table'
include { QC_PEAK_STATS } from './modules/local/qc_peak_stats'
include { QC_TSS_TABLE } from './modules/local/qc_tss_table'
include { TN5_SHIFT } from './modules/local/tn5_shift'
include { TRIM } from './modules/local/trim'
include { TSS_ENRICHMENT } from './modules/local/tss_enrichment'
include { TXDB_BUILD } from './modules/local/txdb_build'

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
            def repl = file(replacement.toString())
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

// Trimmed, non-empty items of a comma-separated value (null-safe).
def csvList(v) {
    if (v == null) return []
    return v.toString().split(',').collect { s -> s.trim() }.findAll { s -> s }
}

// Trimmed text of a config or CSV value ('' for null). Unlike ?:, keeps a numeric 0.
def textOf(v) {
    return v == null ? '' : v.toString().trim()
}

// Read pairs per ALIGN task; 0 = no chunking (blank/null/0). Validated in prepareRun().
def chunkPairs() {
    def t = textOf(params.align_chunk_pairs)
    return t ? t.toLong() : 0L
}

// Number of ALIGN chunks of a library with `pairs` read pairs, `per` pairs each (the last one may be shorter).
def chunkCount(long pairs, long per) { return (pairs + per - 1).intdiv(per).toInteger() }

// Reads a small CSV: UTF-8 BOM and whitespace stripped from keys and values, fully blank rows
// dropped, every required column checked (values may be empty). Fails (exit 1) on a missing file
// or missing columns. Returns a list of row maps.
def readTable(String label, path, List required) {
    def f = file(path.toString())
    if (!f.exists()) fail("${label} not found: ${path}")
    def headerLine = f.readLines().find { l -> l.replace('﻿', '').trim() }
    def header = headerLine ? headerLine.split(',', -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() } : []
    def absent = required.findAll { c -> !(c in header) }
    if (absent) fail("${label} ${path} is missing column(s): ${absent.join(', ')}. Required columns: ${required.join(',')}")
    return f.splitCsv(header: true, quote: '"')
        .collect { row -> row.collectEntries { k, v -> [(k.toString().replace('﻿', '').trim()), (v ?: '').toString().trim()] } }
        .findAll { row -> row.values().any { v -> v } }
}

// Column names of a CSV header (BOM and quotes removed; [] for an empty file).
def csvHeaderOf(path) {
    def line = file(path.toString()).readLines().find { l -> l.replace('﻿', '').trim() }
    return line ? line.split(',', -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() } : []
}

// bowtie2 index files of a prefix (<prefix>*.bt2*), sorted; empty when the index is missing.
def bowtieIndexFiles(String prefix) {
    def base = file(prefix).name
    return files("${prefix}*.bt2*").findAll { f -> f.name.startsWith("${base}.".toString()) }.sort { f -> f.name }
}

// Output prefix of a MACS3 peak group, as the MACS3 module names it (used in the banner only):
// {group}.qval{q}[.minlen{l}][.maxgap{g}][.fe{fe}][.broad]; a part is left out when its param is off.
def macs3NameOf(String group) {
    def p = "${group}.qval${params.macs3_qvalue}"
    if (textOf(params.macs3_min_length)) p += ".minlen${params.macs3_min_length}"
    if (textOf(params.macs3_max_gap)) p += ".maxgap${params.macs3_max_gap}"
    if (textOf(params.macs3_fe_cutoff) && textOf(params.macs3_fe_cutoff).toDouble() > 1) p += ".fe${params.macs3_fe_cutoff}"
    if (isTrue(params.macs3_broad)) p += '.broad'
    return p.toString()
}

// ---- Downstream analyses: region sets and labels ----------------------------------------------------

// Region sets that are differential tables (DIFF_ACCESS, or diff_tables), and the sets DIFF_ANNOTATE writes.
def motifDirectRegions() { return ['promoters', 'peaks'] }
def motifSetRegions() { return ['all_sites', 'promoter_sites', 'non-promoter_sites', 'distal_sites'] }

// Label of a motif region set (the region column of the motif summaries; legacy wording).
def motifLabelOf(String region) {
    def labels = [promoters: 'Promoters', peaks: 'Peaks', all_sites: 'All sites', promoter_sites: 'Promoter sites',
                  'non-promoter_sites': 'Non-promoter sites', distal_sites: 'Distal sites']
    return labels[region] ?: region
}

// BSgenome package (chromosome lengths) and annotation package (gene symbols) of the species and genome build.
// params bsgenome / org_db override them.
def speciesPackages() {
    def mouse = textOf(params.species) == 'Mus musculus'
    def build = textOf(params.genome_build)
    return [bsgenome: textOf(params.bsgenome) ?: ((mouse ? 'BSgenome.Mmusculus.UCSC.' : 'BSgenome.Hsapiens.UCSC.') + build),
            org_db: textOf(params.org_db) ?: (mouse ? 'org.Mm.eg.db' : 'org.Hs.eg.db')]
}

// Downstream analyses (differential accessibility, peak annotation, motifs). Parses comparisons.csv, peak_files.csv,
// diff_tables.csv and motif_targets.csv whenever they are set, and checks them against the switches and the samplesheet.
// Errors are appended to errs (reported by the caller). Returns the plan that workflow {} reads:
//   comps: one map per comparison (name, groupCol, blockCol, pooled, fixedDisp, caseVals, ctrlVals, given (regions from
//          diff_tables), computed (region sets DIFF_ACCESS computes), hasPeaks, ids (controls then cases), ctrlIds, caseIds,
//          blockOf [sample id: block], ctrlPeaks, casePeaks (peak groups, in the listed order))
def prepareDownstream(List errs, List rows, List header, Map peakGroups, String idRe) {
    def runDiff = isTrue(params.run_differential)
    def runAnn = isTrue(params.run_annotation)
    def runMot = isTrue(params.run_motifs)
    def boolOk = ['', 'true', 'false', '1', '0', 'yes', 'no']
    def pkg = speciesPackages()
    def gtf = textOf(params.gtf)

    // ---- Fixed parameters ---------------------------------------------------------------------------
    if (!(textOf(params.species) in ['Homo sapiens', 'Mus musculus'])) errs << "species: '${textOf(params.species)}' must be Homo sapiens or Mus musculus"
    if (!(textOf(params.plot_format) in ['svg', 'pdf'])) errs << "plot_format: '${textOf(params.plot_format)}' must be svg or pdf"
    def diffRegions = csvList(params.diff_regions)
    if (!(new ArrayList(diffRegions).sort() in [['peaks'], ['peaks', 'promoters']])) {
        errs << "diff_regions: '${textOf(params.diff_regions)}' must be peaks or peaks,promoters"
    }
    def dedupText = textOf(params.diff_dedup)
    if (!(dedupText in ['', 'true', 'false', 'auto'])) errs << "diff_dedup: '${dedupText}' must be null, true, false or auto"
    if (!(textOf(params.diff_pe) in ['both', 'none'])) errs << "diff_pe: '${textOf(params.diff_pe)}' must be both or none"
    if (!(textOf(params.diff_table) in ['best', 'combined'])) errs << "diff_table: '${textOf(params.diff_table)}' must be best or combined"
    ['diff_minq', 'diff_max_frag', 'diff_bin_size', 'diff_lfc', 'diff_fdr', 'diff_qlf_min_libs', 'diff_fixed_dispersion',
     'diff_promoter_up', 'diff_promoter_down', 'diff_promoter_min_fc', 'diff_cpus', 'motif_lfc', 'motif_min_score',
     'motif_max_kmer', 'motif_seed', 'motif_top', 'motif_padj_primary', 'motif_padj_relaxed', 'motif_min_bin_peaks',
     'motif_target_min_peaks', 'motif_cpus', 'motif_common_padj', 'motif_common_min_lines'].each { k ->
        if (!isNum(params[k])) errs << "${k}: '${textOf(params[k])}' must be a number"
    }
    if (textOf(params.annotation_fdr) && !isNum(params.annotation_fdr)) errs << "annotation_fdr: '${textOf(params.annotation_fdr)}' must be a number (blank = all regions)"
    if (!(textOf(params.annotation_tss_region) ==~ /-?[0-9]+\s*,\s*-?[0-9]+/)) {
        errs << "annotation_tss_region: '${textOf(params.annotation_tss_region)}' must be two integers, e.g. -2000,200"
    }
    if (!(textOf(params.motif_match) in ['substring', 'exact'])) errs << "motif_match: '${textOf(params.motif_match)}' must be substring or exact"
    def motifRegions = csvList(params.motif_regions)
    if (!motifRegions) errs << "motif_regions is empty"
    motifRegions.findAll { r -> !(r in motifDirectRegions() || r in motifSetRegions()) }.each { r ->
        errs << "motif_regions: '${r}' is not a region set (promoters, peaks, all_sites, promoter_sites, non-promoter_sites, distal_sites)"
    }
    if (textOf(params.motif_pwm_rds) && !file(textOf(params.motif_pwm_rds)).exists()) errs << "motif_pwm_rds not found: ${textOf(params.motif_pwm_rds)}"

    // ---- peak_files.csv: existing peak files, usable as peak groups ------------------------------------
    def peakFiles = []          // [name, path]
    def peakFileNames = []
    if (textOf(params.peak_files)) {
        readTable('peak_files', params.peak_files.toString(), ['name', 'path']).each { r ->
            if (!(r.name ==~ idRe)) errs << "peak_files name '${r.name}' must match ^${idRe}\$ (it is used as a peak group name)"
            else if (r.name in peakFileNames) errs << "duplicate peak_files name: ${r.name}"
            else if (peakGroups.containsKey(r.name)) errs << "peak_files name '${r.name}' is also a peak_group of the samplesheet"
            if (!r.path) errs << "peak_files ${r.name}: path is empty"
            else if (!file(r.path).exists()) errs << "peak_files ${r.name}: file not found: ${r.path}"
            else if (!(r.path ==~ /(?i).*\.(narrowPeak|broadPeak|bed)/)) errs << "peak_files ${r.name}: ${r.path} must end with .narrowPeak, .broadPeak or .bed"
            peakFileNames << r.name
            peakFiles << [name: r.name, path: r.path]
        }
    }

    // ---- comparisons.csv ---------------------------------------------------------------------------------
    def comps = []              // one map per row (see the header of this function)
    def compNames = []
    def compFile = textOf(params.comparisons)
    if (!compFile) {
        if (runDiff || runAnn || runMot) {
            errs << "comparisons is required when run_differential, run_annotation or run_motifs is on (CSV: name,case,control,group_column,block_column,case_peaks,control_peaks,pooled,fixed_dispersion)"
        }
    } else {
        def compCols = ['name', 'case', 'control', 'group_column', 'block_column', 'case_peaks', 'control_peaks', 'pooled', 'fixed_dispersion']
        def unknownC = csvHeaderOf(compFile).findAll { c -> c && !(c in compCols) }
        if (unknownC) errs << "comparisons ${compFile} has unknown column(s): ${unknownC.join(', ')}. Allowed columns: ${compCols.join(', ')}"
        readTable('comparisons', compFile, ['name']).each { r ->
            def name = textOf(r.name)
            if (!(name ==~ idRe)) errs << "comparisons: name '${name}' must match ^${idRe}\$ (it is used in file names)"
            else if (name in compNames) errs << "duplicate comparisons name: ${name}"
            compNames << name
            def gc = textOf(r.group_column) ?: 'condition'
            if (!(gc in header)) errs << "comparisons ${name}: group_column '${gc}' is not a samplesheet column"
            def bc = textOf(r.block_column)
            if (bc && !(bc in header)) errs << "comparisons ${name}: block_column '${bc}' is not a samplesheet column"
            def pooledText = textOf(r.pooled).toLowerCase()
            if (!(pooledText in boolOk)) errs << "comparisons ${name}: pooled '${textOf(r.pooled)}' must be true or false"
            def fd = textOf(r.fixed_dispersion)
            if (fd && !isNum(fd)) errs << "comparisons ${name}: fixed_dispersion '${fd}' must be a number"
            comps << [name: name, groupCol: gc, blockCol: bc, pooled: isTrue(pooledText), fixedDisp: fd ? fd.toDouble() : null,
                      caseVals: csvList(r.case), ctrlVals: csvList(r.control),
                      casePeaksSpec: csvList(r.case_peaks), ctrlPeaksSpec: csvList(r.control_peaks)]
        }
    }

    // ---- diff_tables.csv: differential tables that are given, not computed -------------------------------
    def diffTables = []         // [name, region, path]
    def givenOf = [:]           // comparison name -> regions given in diff_tables
    if (textOf(params.diff_tables)) {
        readTable('diff_tables', params.diff_tables.toString(), ['name', 'regions', 'path']).each { r ->
            if (!(r.name in compNames)) errs << "diff_tables ${r.name}: no comparison of that name in comparisons"
            if (!(r.regions in ['peaks', 'promoters'])) errs << "diff_tables ${r.name}: regions '${r.regions}' must be peaks or promoters"
            def seen = givenOf[r.name] ?: []
            if (r.regions in seen) errs << "diff_tables: ${r.name} / ${r.regions} is listed twice"
            if (!r.path) errs << "diff_tables ${r.name} / ${r.regions}: path is empty"
            else if (!file(r.path).exists()) errs << "diff_tables ${r.name} / ${r.regions}: file not found: ${r.path}"
            givenOf[r.name] = seen + [r.regions]
            diffTables << [name: r.name, region: r.regions, path: r.path]
        }
    }

    // ---- Region sets each comparison computes: diff_regions minus the tables given (only with run_differential) -------
    comps.each { c ->
        c.given = givenOf[c.name] ?: []
        c.computed = runDiff ? diffRegions.findAll { rg -> !(rg in c.given) } : []
        c.hasPeaks = ('peaks' in c.computed) || ('peaks' in c.given)
        c.ids = []
        c.ctrlIds = []
        c.caseIds = []
        c.blockOf = [:]
        c.ctrlPeaks = []
        c.casePeaks = []
    }

    // ---- Samples and peak groups of each comparison --------------------------------------------------------
    comps.each { c ->
        c.caseVals.each { v ->
            if (!rows.any { row -> row[c.groupCol] == v }) errs << "comparisons ${c.name}: case value '${v}' is not a ${c.groupCol} of any active sample"
        }
        c.ctrlVals.each { v ->
            if (!rows.any { row -> row[c.groupCol] == v }) errs << "comparisons ${c.name}: control value '${v}' is not a ${c.groupCol} of any active sample"
        }
        c.caseIds = rows.findAll { row -> row[c.groupCol] in c.caseVals }.collect { row -> row.sample_id }
        c.ctrlIds = rows.findAll { row -> row[c.groupCol] in c.ctrlVals }.collect { row -> row.sample_id }
        def both = c.caseIds.findAll { id -> id in c.ctrlIds }
        if (both) errs << "comparisons ${c.name}: sample(s) on both sides: ${both.join(', ')}"
        if (c.computed) {
            if (!c.caseVals) errs << "comparisons ${c.name}: case is empty (its ${c.computed.join(' and ')} table is computed from case and control)"
            if (!c.ctrlVals) errs << "comparisons ${c.name}: control is empty (its ${c.computed.join(' and ')} table is computed from case and control)"
            if (!c.caseIds) errs << "comparisons ${c.name}: the case side has no active sample (${c.groupCol} in '${c.caseVals.join(',')}')"
            if (!c.ctrlIds) errs << "comparisons ${c.name}: the control side has no active sample (${c.groupCol} in '${c.ctrlVals.join(',')}')"
            if (c.blockCol) {
                (c.ctrlIds + c.caseIds).each { id ->
                    def v = rows.find { row -> row.sample_id == id }[c.blockCol]
                    if (!v) errs << "comparisons ${c.name}: block_column '${c.blockCol}' is empty for sample ${id}"
                    else if (!(v ==~ idRe)) errs << "comparisons ${c.name}: block_column value '${v}' (sample ${id}) must match ^${idRe}\$"
                    else c.blockOf[id] = v
                }
            }
            c.casePeaks = c.casePeaksSpec ?: c.caseIds.collect { id -> rows.find { row -> row.sample_id == id }.peak_group }.findAll { g -> g }.unique()
            c.ctrlPeaks = c.ctrlPeaksSpec ?: c.ctrlIds.collect { id -> rows.find { row -> row.sample_id == id }.peak_group }.findAll { g -> g }.unique()
            if (!c.casePeaks) errs << "comparisons ${c.name}: case_peaks is empty and no case sample has a peak_group"
            if (!c.ctrlPeaks) errs << "comparisons ${c.name}: control_peaks is empty and no control sample has a peak_group"
            (c.casePeaks + c.ctrlPeaks).unique().each { g ->
                if (peakGroups.containsKey(g)) {
                    if (!isTrue(params.run_macs3)) errs << "comparisons ${c.name}: peak group '${g}' is a MACS3 peak group, which needs run_macs3"
                } else if (!(g in peakFileNames)) {
                    errs << "comparisons ${c.name}: peak group '${g}' is neither a peak_group of the samplesheet nor a name in peak_files"
                }
            }
        }
        c.ids = c.ctrlIds + c.caseIds
    }

    // ---- Annotation, TxDb and promoters ----------------------------------------------------------------------
    if (runAnn) {
        if (!gtf) errs << "run_annotation needs gtf (the annotation TxDb is built from it)"
        comps.findAll { c -> !c.hasPeaks }.each { c ->
            errs << "run_annotation: comparison ${c.name} has no peaks table (run_differential computes it, or give it in diff_tables)"
        }
    }
    def promotersComputed = comps.any { c -> 'promoters' in c.computed }
    def needTxdb = runAnn || promotersComputed
    if (needTxdb && !gtf) errs << "gtf is required when promoters are computed or run_annotation is on (the TxDb is built from it)"

    // ---- Motifs: every comparison needs the table of each motif_regions set ----------------------------------------
    if (runMot) {
        motifRegions.each { r ->
            comps.each { c ->
                def ok = (r in motifDirectRegions()) ? ((r in c.computed) || (r in c.given)) : runAnn
                if (!ok) {
                    def how = (r in motifDirectRegions()) ? "run_differential with ${r} in diff_regions, or ${r} in diff_tables" : 'run_annotation'
                    errs << "run_motifs: comparison ${c.name} has no ${r} table (needed by motif_regions; ${how})"
                }
            }
        }
    }

    // ---- motif_targets.csv: target-vs-other analyses on the all_sites table of their comparison ----------------------
    def motifTargets = []       // [name, type, target, comparison, label, minPeaks, bed]
    if (textOf(params.motif_targets)) {
        if (!runAnn) errs << "motif_targets needs run_annotation (a target analysis runs on the all_sites table of its comparison)"
        def tNames = []
        readTable('motif_targets', textOf(params.motif_targets), ['name', 'type', 'target', 'comparison']).each { r ->
            def name = textOf(r.name)
            if (!(name ==~ idRe)) errs << "motif_targets name '${name}' must match ^${idRe}\$ (it is used in file names)"
            else if (name in tNames) errs << "duplicate motif_targets name: ${name}"
            else if (name in compNames) errs << "motif_targets ${name}: same name as a comparison (their motif folders would clash)"
            tNames << name
            if (!(r.type in ['gene', 'regions'])) errs << "motif_targets ${name}: type '${r.type}' must be gene or regions"
            if (!r.target) errs << "motif_targets ${name}: target is empty"
            else if (r.type == 'gene' && !(r.target ==~ /[A-Za-z0-9_.\-]+/)) errs << "motif_targets ${name}: gene symbol '${r.target}' may contain only letters, digits, _ . -"
            else if (r.type == 'regions' && !file(r.target).exists()) errs << "motif_targets ${name}: BED not found: ${r.target}"
            if (!(r.comparison in compNames)) errs << "motif_targets ${name}: comparison '${r.comparison}' is not in comparisons"
            def mp = textOf(r.min_peaks)
            if (mp && !(mp ==~ /[0-9]+/)) errs << "motif_targets ${name}: min_peaks '${mp}' must be a whole number"
            def label = textOf(r.label) ?: (r.type == 'gene' ? "${r.target} peaks" : 'Target regions')
            if ((label =~ /["`\\$]/).find()) errs << "motif_targets ${name}: label '${label}' may not contain quotes, backticks, backslashes or dollar signs"
            motifTargets << [name: name, type: r.type, target: r.target, comparison: r.comparison, label: label,
                             minPeaks: mp ?: textOf(params.motif_target_min_peaks),
                             bed: (r.type == 'regions' && r.target) ? file(r.target) : null]
        }
    }

    return [runDiff: runDiff, runAnn: runAnn, runMot: runMot, needTxdb: needTxdb, comps: comps, compNames: compNames,
            pooledNames: comps.findAll { c -> c.pooled }.collect { c -> c.name },
            peakFiles: peakFiles, diffTables: diffTables, motifTargets: motifTargets, motifRegions: motifRegions,
            diffRegions: diffRegions, diffDedup: dedupText ?: (isTrue(params.remove_duplicates) ? 'false' : 'true'),
            bsgenome: pkg.bsgenome, orgDb: pkg.org_db]
}

def prepareRun() {
    if (params.help) {
        log.info """
        ATACSEQ pipeline (ATAC-seq, paired-end)
        ---------------------------------------
        Usage:
          nextflow run ${workflow.projectDir} -profile pmacs -params-file config.yaml

        Required params (set in config.yaml):
          samplesheet     CSV: sample_id,fastq_1,fastq_2,condition,replicate,tmm_group,tmm_group_legacy,peak_group,exclude,bam
                          (sample_id required; each active row has fastq_1 and fastq_2, or a final BAM in bam, not both;
                          the other columns are optional, blank when absent)
          genome_index    bowtie2 index prefix of the genome

        Optional params:
          tmm_group_column       samplesheet column that defines the TMM groups (default tmm_group; or tmm_group_legacy)
          peak_merges            CSV: name,peak_groups (comma list, quoted): one merged BED per row (needs run_macs3)
          blacklist              blacklist BED (empty = skip BLACKLIST_FILTER)
          contamination_index    bowtie2 index prefix of the contaminant genomes (empty = skip CONTAMINATION)
          tss_bed                TSS BED (needed when run_tss_enrichment is on and some sample has a tmm_group)
          effective_genome_size  MACS3 -g (needed when run_macs3 is on and some sample has a peak_group)
          chrom_sizes, gtf       recorded and validated only (gtf is also needed by promoters and by run_annotation)
          align_chunk_pairs      read pairs per ALIGN task (default 10000000; 0 or blank = one task per library)
          run_fastqc, trim, run_tss_enrichment, run_macs3, run_multiqc, tn5_shift, remove_duplicates
          project_name, outdir (./results)

        Downstream analyses (all off by default):
          run_differential       DIFF_ACCESS per comparison (csaw/edgeR; needs comparisons)
          run_annotation         DIFF_ANNOTATE per comparison (ChIPseeker; needs gtf and a peaks table)
          run_motifs             MOTIF_ENRICHMENT per comparison and motif_regions set and per motif_targets row; MOTIF_COMMON
          comparisons            CSV: name,case,control,group_column,block_column,case_peaks,control_peaks,pooled,fixed_dispersion
          peak_files             CSV: name,path (peak files usable as peak groups)
          diff_tables            CSV: name,regions,path (differential tables given, not computed)
          motif_targets          CSV: name,type,target,comparison,label,min_peaks (gene | regions target vs other)
          diff_*, annotation_*, motif_*, species, genome_build, bsgenome, org_db, plot_font, plot_format: see nextflow.config

        The other keys and their defaults are in nextflow.config. Replace one step's script with
        <script name>: path in config.yaml (see stepScripts()).
        """
        System.exit(0)
    }

    def errs = []
    def idRe = '[A-Za-z0-9][A-Za-z0-9_.-]*'
    def boolOk = ['', 'true', 'false', '1', '0', 'yes', 'no']

    // ---- Top-level params ---------------------------------------------------
    def allowed = [bowtie2_mode: ['local', 'end_to_end'], bigwig_units: ['none', 'cpm', 'rpkm']]
    allowed.each { k, vals ->
        def v = textOf(params[k])
        if (!(v in vals)) errs << "${k}: '${v}' is not valid (one of: ${vals.join(', ')})"
    }
    // Numbers: required ones must be numeric; the null-ok ones may be blank (= off).
    ['trim_minlen', 'bowtie2_maxins', 'alignment_min_pct', 'tss_flank', 'tmm_bin_size', 'bigwig_bin_size',
     'macs3_qvalue', 'macs3_fe_cutoff', 'macs3_broad_cutoff'].each { k ->
        if (!isNum(params[k])) errs << "${k}: '${textOf(params[k])}' must be a number"
    }
    ['bowtie2_minins', 'macs3_min_length', 'macs3_max_gap'].each { k ->
        if (textOf(params[k]) != '' && !isNum(params[k])) errs << "${k}: '${textOf(params[k])}' must be a number (blank = off)"
    }
    def chunkText = textOf(params.align_chunk_pairs)
    if (chunkText != '' && !(chunkText ==~ /[0-9]+/)) errs << "align_chunk_pairs: '${chunkText}' must be an integer >= 0 (0 or blank = no chunking)"
    if (!isNum(params.mapq)) errs << "mapq: '${textOf(params.mapq)}' must be a number >= 0"
    else if (params.mapq.toString().trim().toDouble() < 0) errs << "mapq: '${textOf(params.mapq)}' must be a number >= 0"
    def tmmMax = textOf(params.tmm_max_frag)
    if (!(tmmMax == 'auto' || (tmmMax ==~ /[0-9]+/ && tmmMax.toInteger() > 0))) {
        errs << "tmm_max_frag: '${tmmMax}' must be 'auto' or a positive integer"
    }
    ['tn5_shift', 'remove_duplicates', 'trim', 'run_fastqc', 'run_macs3', 'run_tss_enrichment', 'run_multiqc',
     'macs3_call_summits', 'macs3_broad', 'bigwig_extend_reads', 'bigwig_center_reads',
     'run_differential', 'run_annotation', 'run_motifs'].each { k ->
        if (!(textOf(params[k]).toLowerCase() in boolOk)) errs << "${k}: '${textOf(params[k])}' must be true or false"
    }
    if (!(textOf(params.lsf_exclude_hosts) ==~ /[A-Za-z0-9_.\-, ]*/))
        errs << "lsf_exclude_hosts: '${textOf(params.lsf_exclude_hosts)}' must be a comma list of host names"
    failList("invalid config.yaml value(s)", errs)

    if (!textOf(params.samplesheet)) fail("missing required param: samplesheet (set in your config.yaml).")

    // ---- Samplesheet ----------------------------------------------------------
    def sheet = params.samplesheet.toString()
    def sheetCols = ['sample_id', 'fastq_1', 'fastq_2', 'condition', 'replicate', 'tmm_group', 'tmm_group_legacy', 'peak_group', 'exclude', 'bam']
    def allRows = readTable('samplesheet', sheet, ['sample_id', 'fastq_1', 'fastq_2'])
    if (!allRows) fail("samplesheet has no data rows: ${sheet}")
    def headerLine = file(sheet).readLines().find { l -> l.replace('﻿', '').trim() }
    def header = headerLine.split(',', -1).collect { h -> h.replace('﻿', '').replace('"', '').trim() }
    def unknown = header.findAll { c -> c && !(c in sheetCols) }
    if (unknown) errs << "samplesheet ${sheet} has unknown column(s): ${unknown.join(', ')}. Allowed columns: ${sheetCols.join(', ')}"
    // An absent optional column is blank for every row.
    allRows.each { row -> sheetCols.each { c -> row[c] = row[c] ?: '' } }

    def tmmCol = textOf(params.tmm_group_column)
    if (!(tmmCol in header)) errs << "tmm_group_column: '${tmmCol}' is not a column of ${sheet} (columns: ${header.join(', ')})"

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

    // Only the rows that are not excluded are processed; their group values define the groups.
    def rows = allRows.findAll { row -> !row.excluded }
    def excludedIds = allRows.findAll { row -> row.excluded }.collect { row -> row.sample_id }
    if (!rows) fail("every row is excluded: ${sheet}")
    def ids = rows.collect { row -> row.sample_id }
    // Libraries with FASTQ (the FASTQ steps) and libraries with a final BAM (BAM_PREALIGNED); a row has one or the other.
    def fastqRows = rows.findAll { row -> !row.bam }
    def bamRows = rows.findAll { row -> row.bam }
    def exGroup = allRows.findAll { row -> row.excluded && (row[tmmCol] || row.peak_group) }.collect { row -> row.sample_id }
    if (exGroup) log.info "excluded rows with a group value (ignored): ${exGroup.join(', ')}"

    // Input files and FASTQ names.
    rows.each { row ->
        if (row.bam) {
            if (row.fastq_1 || row.fastq_2) errs << "${row.sample_id}: bam is set, so fastq_1 and fastq_2 must be empty (a library has FASTQ or a BAM, not both)"
            if (!row.bam.endsWith('.bam')) errs << "${row.sample_id}: bam '${row.bam}' must be a file ending in .bam"
            else if (!file(row.bam).exists()) errs << "${row.sample_id}: bam not found: ${row.bam}"
        } else {
            ['fastq_1', 'fastq_2'].each { k ->
                if (!row[k]) errs << "${row.sample_id}: ${k} is empty"
                else if (!file(row[k]).exists()) errs << "${row.sample_id}: ${k} not found: ${row[k]}"
            }
            if (row.fastq_1 && row.fastq_1 == row.fastq_2) errs << "${row.sample_id}: fastq_1 and fastq_2 are the same file (${row.fastq_1})"
        }
    }
    def fqPaths = fastqRows.collectMany { row -> [row.fastq_1, row.fastq_2] }.findAll { p -> p }
    def fqNames = fqPaths.collect { p -> file(p).name }
    fqNames.findAll { n -> fqNames.count(n) > 1 }.unique().each { n ->
        errs << "two FASTQ files share the name ${n}: " + fqPaths.findAll { p -> file(p).name == n }.join(', ')
    }

    // Group values: TMM column and peak_group must be usable in file names.
    rows.each { row ->
        [tmmCol, 'peak_group'].unique().each { col ->
            def v = row[col]
            if (v && !(v ==~ idRe)) errs << "${row.sample_id}: ${col} '${v}' must match ^${idRe}\$ (it is used in file names)"
        }
    }
    def tmmOf = [:]
    def tmmGroups = [:]     // group -> [sample ids], samplesheet order
    def peakOf = [:]
    def peakGroups = [:]    // group -> [sample ids], samplesheet order
    rows.each { row ->
        def t = row[tmmCol]
        if (t) {
            tmmOf[row.sample_id] = t
            if (!tmmGroups.containsKey(t)) tmmGroups[t] = []
            tmmGroups[t] << row.sample_id
        }
        def p = row.peak_group
        if (p) {
            peakOf[row.sample_id] = p
            if (!peakGroups.containsKey(p)) peakGroups[p] = []
            peakGroups[p] << row.sample_id
        }
    }
    tmmGroups.findAll { g, sids -> sids.size() == 1 }.each { g, sids ->
        errs << "tmm group '${g}' has only one sample (${sids[0]}); TMM needs at least 2"
    }
    failList("${sheet}: invalid sample / group values", errs)

    // Paths are resolved only now: every file in the sheet exists.
    def fqById = [:]
    fastqRows.each { row -> fqById[row.sample_id] = [file(row.fastq_1), file(row.fastq_2)] }

    // ---- peak_merges ----------------------------------------------------------
    def merges = []   // [name, [peak groups]]
    if (textOf(params.peak_merges)) {
        def mrows = readTable('peak_merges', params.peak_merges.toString(), ['name', 'peak_groups'])
        def mNames = mrows.collect { r -> r.name }
        def dupM = mNames.findAll { n -> mNames.count(n) > 1 }.unique()
        if (dupM) errs << "duplicate peak_merges name(s): ${dupM.join(', ')}"
        mrows.each { r ->
            if (!(r.name ==~ idRe)) errs << "peak_merges name '${r.name}' must match ^${idRe}\$ (it is used in file names)"
            def gs = csvList(r.peak_groups)
            if (!gs) errs << "peak_merges ${r.name}: peak_groups is empty"
            gs.findAll { g -> !peakGroups.containsKey(g) }.each { g ->
                errs << "peak_merges ${r.name}: peak group '${g}' has no active sample in ${sheet}"
            }
            merges << [r.name, gs]
        }
        if (!isTrue(params.run_macs3)) errs << "peak_merges is set but run_macs3 is false"
    }
    failList("peak_merges: invalid", errs)

    // ---- Reference files ------------------------------------------------------
    def refErrs = []
    def prefix = textOf(params.genome_index)
    def indexFiles = []
    def indexBase = null
    if (!prefix) {
        refErrs << "genome_index is required (bowtie2 index prefix)"
    } else {
        indexBase = file(prefix).name
        indexFiles = bowtieIndexFiles(prefix)
        if (!indexFiles) refErrs << "genome_index: no bowtie2 index files (${prefix}.1.bt2 or ${prefix}.1.bt2l ...) found for prefix ${prefix}"
    }
    def blacklist = textOf(params.blacklist)
    if (blacklist && !file(blacklist).exists()) refErrs << "blacklist not found: ${blacklist}"
    ['chrom_sizes', 'gtf'].each { k ->
        def v = textOf(params[k])
        if (v && !file(v).exists()) refErrs << "${k} not found: ${v}"
    }
    def tss = textOf(params.tss_bed)
    if (isTrue(params.run_tss_enrichment) && tmmGroups) {
        if (!tss) refErrs << "tss_bed is required (run_tss_enrichment is on)"
        else if (!file(tss).exists()) refErrs << "tss_bed not found: ${tss}"
    }
    def egs = textOf(params.effective_genome_size)
    if (isTrue(params.run_macs3)) {
        if (!egs) {
            if (peakGroups) refErrs << "effective_genome_size is required (run_macs3 is on)"
        } else if (!isNum(egs)) {
            refErrs << "effective_genome_size: '${egs}' must be a number"
        }
    }
    def contamPrefix = textOf(params.contamination_index)
    def contamFiles = []
    def contamBase = null
    if (contamPrefix) {
        contamBase = file(contamPrefix).name
        contamFiles = bowtieIndexFiles(contamPrefix)
        if (!contamFiles) refErrs << "contamination_index: no bowtie2 index files found for prefix ${contamPrefix}"
    }
    failList("reference file problems", refErrs)
    def refs = [index_files: indexFiles, index_base: indexBase,
                contam_files: contamFiles, contam_base: contamBase,
                blacklist: blacklist ? file(blacklist) : null,
                chrom_sizes: textOf(params.chrom_sizes) ? file(textOf(params.chrom_sizes)) : null,
                gtf: textOf(params.gtf) ? file(textOf(params.gtf)) : null,
                tss_bed: tss ? file(tss) : null]

    // ---- Downstream analyses (their own checks; see prepareDownstream) -----------------
    def downstream = prepareDownstream(errs, rows, header, peakGroups, idRe)
    failList("downstream analyses: invalid value(s)", errs)

    // ---- Summary --------------------------------------------------------------
    def steps = []
    ['run_fastqc', 'trim', 'remove_duplicates', 'tn5_shift', 'run_tss_enrichment', 'run_macs3', 'run_multiqc',
     'run_differential', 'run_annotation', 'run_motifs'].each { k ->
        if (isTrue(params[k])) steps << k
    }
    if (refs.blacklist != null) steps << 'blacklist_filter'
    if (refs.contam_files) steps << 'contamination'
    if (textOf(params.remove_chroms)) steps << "remove_chroms (${textOf(params.remove_chroms)})"
    if (merges) steps << 'peak_merges'
    if (bamRows) steps << "pre-aligned BAM (${bamRows.size()})"
    def tmmText = tmmGroups.collect { g, sids -> "${g} (${sids.size()})" }.join(', ')
    def peakText = peakGroups.collect { g, sids -> "${g} (${sids.size()})" }.join(', ')
    def mergeText = merges.collect { m -> "${m[0]} = ${m[1].join('+')}" }.join(', ')
    def alignText = chunkPairs() > 0 ? "bowtie2 in chunks of ${chunkPairs()} pairs" : 'bowtie2, one task per library'
    def downText = (downstream.runDiff || downstream.runAnn || downstream.runMot)
        ? "comparisons ${downstream.compNames.size()} (${downstream.compNames.join(', ')}); diff regions ${downstream.diffRegions.join(',')}; motif regions ${downstream.motifRegions.join(',')}; motif targets ${downstream.motifTargets.size()}"
        : 'off'
    log.info """
    ATACSEQ  ~  ${params.project_name}
    -----------------------------------------
      libraries      : ${ids.size()} active (${fastqRows.size()} FASTQ, ${bamRows.size()} final BAM), ${excludedIds.size()} excluded${excludedIds ? ' (' + excludedIds.join(', ') + ')' : ''}
      alignment      : ${alignText}
      TMM groups     : ${tmmGroups ? tmmText : 'none'} (column: ${tmmCol})
      peak groups    : ${peakGroups ? peakText : 'none'}
      peak merges    : ${merges ? mergeText : 'none'}
      MACS3 prefix   : ${peakGroups ? 'e.g. ' + macs3NameOf(peakGroups.keySet().toList()[0]) : '(no peak groups)'}
      downstream     : ${downText}
      steps on       : ${steps.join(', ')}
      outdir         : ${params.outdir}
    """

    return [scripts: stepScripts(), images: processContainers(), refs: refs, rows: rows, excludedIds: excludedIds,
            ids: ids, fastqIds: fastqRows.collect { row -> row.sample_id }, bamRows: bamRows, fqById: fqById,
            tmmGroups: tmmGroups, tmmOf: tmmOf, peakGroups: peakGroups, peakOf: peakOf,
            merges: merges, downstream: downstream]
}

// Channel of tuple(group, ids, bams, bais), one per group, emitted as soon as every sample of that group has an analysis BAM.
// members = [group: [sample ids, samplesheet order]], groupOf = [sample id: group]. ids/bams/bais come out in samplesheet
// order, so a task gets the same inputs (and writes the same outputs) as when all BAMs were collected first. Every member of a
// group is an active samplesheet row (prepareRun builds the groups from those rows only), so each one reaches the channel.
def groupedBams(ch, Map groupOf, Map members) {
    return ch.filter { id, bam, bai -> groupOf.containsKey(id) }
        .map { id, bam, bai -> tuple(groupKey(groupOf[id], members[groupOf[id]].size()), id, bam, bai) }
        .groupTuple()
        .map { key, ids, bams, bais ->
            def g = key.toString()
            def pos = members[g].collect { sid -> ids.indexOf(sid) }
            tuple(g, members[g], pos.collect { i -> bams[i] }, pos.collect { i -> bais[i] })
        }
}

workflow {

    def cfg = prepareRun()

    // Each step's script: bin/<name>, or this project's replacement (see stepScripts()).
    def sc = cfg.scripts
    def samplesheet = Channel.value(file(params.samplesheet))

    // Channels that optional blocks add to (all stay empty when their block is off or has no FASTQ library).
    def fastqc_zips = Channel.empty()
    def trim_logs = Channel.empty()
    def macs_xls = Channel.empty()
    def align_logs = Channel.empty()        // S.bowtie2.log (FASTQ libraries)
    def dedup_logs = Channel.empty()        // DEDUP.out.log
    def raw_flagstat = Channel.empty()      // BAM_FILTER.out.raw_flagstat
    def raw_idxstats = Channel.empty()      // BAM_FILTER.out.raw_idxstats
    // Final BAMs of every library: BAM_FINALIZE (FASTQ libraries), mixed with BAM_PREALIGNED (bam rows).
    def final_bams = Channel.empty()        // tuple(sample_id, S.final.bam, S.final.bam.bai)
    def final_flagstat = Channel.empty()    // S.final.flagstat.txt
    def final_idxstats = Channel.empty()    // S.final.idxstats.txt
    def macs_peaks = Channel.empty()        // tuple(peak group, MACS3 peak file), when MACS3 runs

    // ---- A. Preprocess: FASTQ libraries only (skipped when every library is a pre-aligned BAM) ----------------
    if (cfg.fastqIds) {
        def reads = Channel.fromList(cfg.fastqIds.collect { id -> tuple(id, cfg.fqById[id][0], cfg.fqById[id][1]) })
        if (isTrue(params.run_fastqc)) {
            FASTQC_RAW(Channel.fromList(cfg.fastqIds.collectMany { id -> [tuple(id, cfg.fqById[id][0]), tuple(id, cfg.fqById[id][1])] }), 'raw', sc.fastqc)
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
        // ALIGN: with align_chunk_pairs > 0 each library is counted (COUNT_PAIRS), split into chunks of that many pairs
        // (bowtie2 -s/-u, no temporary FASTQs), aligned in parallel and merged (ALIGN_MERGE) as soon as ITS chunks are done.
        def chunk_pairs = chunkPairs()
        def align_bams = null   // tuple(sample_id, S.aligned.bam)
        if (chunk_pairs > 0) {
            COUNT_PAIRS(reads, sc.count_pairs)
            def pair_counts = COUNT_PAIRS.out.count.map { id, f -> tuple(id, f.text.trim().toLong()) }
            def chunk_in = reads.join(pair_counts, failOnMismatch: true, failOnDuplicate: true)
                .flatMap { id, r1, r2, n -> (0..<chunkCount(n, chunk_pairs)).collect { k -> tuple(id, k, k * chunk_pairs, chunk_pairs, r1, r2) } }
            ALIGN(chunk_in, Channel.value(cfg.refs.index_files), cfg.refs.index_base, params.bowtie2_mode, textOf(params.bowtie2_minins),
                  params.bowtie2_maxins, textOf(params.bowtie2_args), sc.align)
            def n_chunks = pair_counts.map { id, n -> tuple(id, chunkCount(n, chunk_pairs)) }
            // groupKey(id, n chunks): the merge fires when that library's chunks are all done; lists are sorted by chunk index.
            def merge_in = ALIGN.out.aligned.combine(n_chunks, by: 0)
                .map { id, k, bam, log, nc -> tuple(groupKey(id, nc), k, bam, log) }
                .groupTuple()
                .map { key, ks, bams, logs ->
                    def order = (0..<ks.size()).toList().sort { i -> ks[i] }
                    tuple(key.toString(), order.collect { i -> bams[i] }, order.collect { i -> logs[i] })
                }
            ALIGN_MERGE(merge_in, sc.merge_bowtie2_logs, sc.align_merge)
            align_bams = ALIGN_MERGE.out.bam
            align_logs = ALIGN_MERGE.out.log
        } else {
            ALIGN(reads.map { id, r1, r2 -> tuple(id, -1, 0, 0, r1, r2) }, Channel.value(cfg.refs.index_files), cfg.refs.index_base,
                  params.bowtie2_mode, textOf(params.bowtie2_minins), params.bowtie2_maxins, textOf(params.bowtie2_args), sc.align)
            align_bams = ALIGN.out.aligned.map { id, k, bam, log -> tuple(id, bam) }
            align_logs = ALIGN.out.aligned.map { id, k, bam, log -> log }
        }
        ALIGNMENT_SUMMARY(align_logs.collect(), params.alignment_min_pct, sc.alignment_summary)
        BAM_FILTER(align_bams, params.mapq, textOf(params.remove_chroms), sc.bam_filter)
        raw_flagstat = BAM_FILTER.out.raw_flagstat
        raw_idxstats = BAM_FILTER.out.raw_idxstats
        if (textOf(params.remove_chroms)) {
            QC_CHROM_FRACTION(BAM_FILTER.out.raw_idxstats.collect(), textOf(params.remove_chroms), sc.qc_chrom_fraction)
        }
        DEDUP(BAM_FILTER.out.bam, isTrue(params.remove_duplicates).toString(), sc.dedup)
        dedup_logs = DEDUP.out.log
        def pre_final = DEDUP.out.bam
        if (cfg.refs.blacklist != null) {
            BLACKLIST_FILTER(DEDUP.out.bam, Channel.value(cfg.refs.blacklist), sc.blacklist_filter)
            pre_final = BLACKLIST_FILTER.out.bam
        }
        BAM_FINALIZE(pre_final, sc.bam_finalize)
        final_bams = BAM_FINALIZE.out.bam
        final_flagstat = BAM_FINALIZE.out.flagstat
        final_idxstats = BAM_FINALIZE.out.idxstats
        if (cfg.refs.contam_files) {
            CONTAMINATION(align_bams, Channel.value(cfg.refs.contam_files), cfg.refs.contam_base, sc.contamination)
        }
    }

    // ---- A2. Pre-aligned BAMs (samplesheet column bam): the final-BAM layout, mixed with BAM_FINALIZE's outputs ----
    // A pre-aligned BAM is assumed final (filtered, duplicate-removed, blacklist-filtered): it skips BAM_FILTER to BLACKLIST_FILTER.
    if (cfg.bamRows) {
        BAM_PREALIGNED(Channel.fromList(cfg.bamRows.collect { row -> tuple(row.sample_id, file(row.bam)) }), sc.bam_prealigned)
        final_bams = final_bams.mix(BAM_PREALIGNED.out.bam)
        final_flagstat = final_flagstat.mix(BAM_PREALIGNED.out.flagstat)
        final_idxstats = final_idxstats.mix(BAM_PREALIGNED.out.idxstats)
    }

    // Analysis BAMs: S.tn5.bam when tn5_shift is on, else S.final.bam. QC, TMM, bigWig, MACS3 and FRiP all read these.
    def analysis_bams = final_bams
    if (isTrue(params.tn5_shift)) {
        TN5_SHIFT(final_bams, sc.tn5_shift)
        analysis_bams = TN5_SHIFT.out.bam
    }

    // ---- B. QC: fragment sizes ----------------------------------------------------
    QC_FRAGMENT_SIZES_TABLE(analysis_bams.flatMap { id, bam, bai -> [bam, bai] }.collect(), samplesheet, sc.qc_fragment_sizes_table)
    QC_FRAGMENT_SIZES_PLOT(QC_FRAGMENT_SIZES_TABLE.out.summary, QC_FRAGMENT_SIZES_TABLE.out.distribution,
                           Channel.value(sc.my_theme), sc.qc_fragment_sizes_plot)

    // ---- C. TMM scale factors, bigWigs, TSS enrichment ------------------------------
    // Samples without a tmm_group get no bigWig and no TSS enrichment.
    if (cfg.tmmGroups) {
        def tmm_in = groupedBams(analysis_bams, cfg.tmmOf, cfg.tmmGroups).map { g, ids, bams, bais -> tuple(g, ids.join(','), bams, bais) }
        NORM_TMM(tmm_in, params.tmm_bin_size, params.mapq, params.tmm_max_frag, params.bigwig_units, sc.norm_tmm)

        // tuple(sample_id, bigwig_multiplier, method) from every group's scale-factor table
        def factors = NORM_TMM.out.sf.flatMap { g, tsv ->
            tsv.splitCsv(header: true, sep: '\t').collect { r -> tuple(r.sample_id, r.bigwig_multiplier, r.method) }
        }
        // Keyed join: every sample with a tmm_group has exactly one scale factor.
        def bw_in = analysis_bams.filter { id, bam, bai -> cfg.tmmOf.containsKey(id) }
            .join(factors, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, bam, bai, multiplier, method)
        BIGWIG(bw_in, params.bigwig_bin_size, params.bigwig_units, isTrue(params.bigwig_extend_reads).toString(),
               isTrue(params.bigwig_center_reads).toString(), params.effective_genome_size, sc.bigwig)

        if (isTrue(params.run_tss_enrichment)) {
            TSS_ENRICHMENT(BIGWIG.out.bw, Channel.value(cfg.refs.tss_bed), params.tss_flank, sc.qc_tss_enrichment)
            QC_TSS_TABLE(TSS_ENRICHMENT.out.tsv.collect(), samplesheet, sc.qc_tss_table)
        }
    }

    // ---- D. Peaks: MACS3 per peak group, peak merges, FRiP, peak statistics ---------
    def peak_files = Channel.empty()   // every peak file for QC_PEAK_STATS (MACS3 peaks and merged BEDs)
    if (isTrue(params.run_macs3) && cfg.peakGroups) {
        // Replicates pooled, no control: one task per peak group.
        def macs_in = groupedBams(analysis_bams, cfg.peakOf, cfg.peakGroups).map { g, ids, bams, bais -> tuple(g, bams) }
        MACS3(macs_in, params.effective_genome_size, params.macs3_qvalue, textOf(params.macs3_min_length), textOf(params.macs3_max_gap),
              params.macs3_fe_cutoff, isTrue(params.macs3_call_summits).toString(), isTrue(params.macs3_broad).toString(),
              params.macs3_broad_cutoff, sc.callpeaks_macs3)
        macs_xls = MACS3.out.xls
        macs_peaks = MACS3.out.peaks
        peak_files = MACS3.out.peaks.map { g, f -> f }

        if (cfg.merges) {
            // Each merge waits only for the peaks of its own groups; files come in the merge's group order.
            def merge_in = MACS3.out.peaks
                .flatMap { g, f -> cfg.merges.findAll { m -> g in m[1] }.collect { m -> tuple(groupKey(m[0], m[1].unique().size()), g, f) } }
                .groupTuple()
                .map { key, gs, fs ->
                    def name = key.toString()
                    def groups = cfg.merges.find { m -> m[0] == name }[1].unique()
                    tuple(name, groups.collect { g -> fs[gs.indexOf(g)] })
                }
            PEAK_MERGE(merge_in, sc.peak_merge)
            peak_files = peak_files.mix(PEAK_MERGE.out.bed.map { n, f -> f })
        }

        // FRiP: each sample against ITS peak group's peak file (waits only for that group's MACS3 peaks).
        // The flagstat is keyed by sample id from its file name.
        def flagstat_keyed = final_flagstat.map { f -> tuple(f.name.replace('.final.flagstat.txt', ''), f) }
        def frip_in = analysis_bams.join(flagstat_keyed, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, bam, bai, flagstat)
            .filter { id, bam, bai, fs -> cfg.peakOf.containsKey(id) }
            .map { id, bam, bai, fs -> tuple(cfg.peakOf[id], id, bam, bai, fs) }
            .combine(MACS3.out.peaks, by: 0)
            .map { g, id, bam, bai, fs, pf -> tuple(id, bam, bai, pf, 'macs3', fs) }
        FRIP(frip_in, sc.qc_frip)
        QC_FRIP_TABLE(FRIP.out.tsv.collect(), samplesheet, sc.qc_frip_table)
        QC_PEAK_STATS(peak_files.collect(), sc.qc_peak_stats)
    }

    // ---- E. Downstream: differential accessibility, peak annotation, motifs (all off by default) ----------------
    // Nothing below is instantiated unless its switch is on. Every gathering is keyed (groupKey with the exact count, or a
    // join on an explicit key); the one unkeyed gather is MOTIF_COMMON, explained at its call.
    def dn = cfg.downstream
    def no_file = file("${projectDir}/assets/NO_FILE")   // placeholder for a path input that a task does not use
    def diff_out = Channel.fromList(dn.diffTables.collect { t -> tuple(t.name, t.region, file(t.path)) })   // tuple(name, region, csv)
    def annotation_out = Channel.empty()                                                                 // tuple(name, set, csv)
    if (dn.runDiff || dn.runAnn || dn.runMot) {
        // TxDb (gene and transcript models from the GTF): promoters and annotation. Its inputs are values, so its output is a
        // value channel that every task can take.
        def txdb = Channel.value(no_file)
        if (dn.needTxdb) {
            TXDB_BUILD(file(textOf(params.gtf)), params.species, sc.txdb_build)
            txdb = TXDB_BUILD.out.txdb
        }

        if (dn.runDiff && dn.comps.any { c -> c.computed }) {
            def dcomps = dn.comps.findAll { c -> c.computed }
            // BAMs of each comparison: controls first, then cases, samplesheet order within a side (groupKey = the member count).
            def bam_in = analysis_bams.flatMap { id, bam, bai ->
                    dcomps.findAll { c -> id in c.ids }.collect { c -> tuple(groupKey(c.name, c.ids.size()), id, bam, bai) }
                }
                .groupTuple()
                .map { key, ids, bams, bais ->
                    def c = dcomps.find { x -> x.name == key.toString() }
                    def pos = c.ids.collect { sid -> ids.indexOf(sid) }
                    tuple(c.name, pos.collect { i -> bams[i] }, pos.collect { i -> bais[i] })
                }
            // Peak files of each side: the MACS3 peak files of the named groups and the peak_files entries, in the listed order
            // (groupKey = the number of listed groups of that side). Key 'name:control' or 'name:case'.
            def peak_items = []   // [key, index, group]
            dcomps.each { c ->
                c.ctrlPeaks.eachWithIndex { g, i -> peak_items << [key: "${c.name}:control", idx: i, group: g] }
                c.casePeaks.eachWithIndex { g, i -> peak_items << [key: "${c.name}:case", idx: i, group: g] }
            }
            def peak_sizes = peak_items.countBy { it.key }
            def peak_src = macs_peaks.mix(Channel.fromList(dn.peakFiles.collect { p -> tuple(p.name, file(p.path)) }))
            def peak_in = peak_src.flatMap { g, f ->
                    peak_items.findAll { it.group == g }.collect { it -> tuple(groupKey(it.key, peak_sizes[it.key]), it.idx, f) }
                }
                .groupTuple()
                .map { key, idxs, fs ->
                    def order = (0..<idxs.size()).toList().sort { i -> idxs[i] }
                    tuple(key.toString(), order.collect { i -> fs[i] })
                }
            def ctrl_peaks = peak_in.filter { k, fs -> k.endsWith(':control') }.map { k, fs -> tuple(k.take(k.lastIndexOf(':')), fs) }
            def case_peaks = peak_in.filter { k, fs -> k.endsWith(':case') }.map { k, fs -> tuple(k.take(k.lastIndexOf(':')), fs) }

            // One task per comparison: BAMs, sample ids, groups, blocks, peaks, and the region sets and dispersion of THIS comparison.
            def diff_in = bam_in.join(ctrl_peaks, failOnMismatch: true, failOnDuplicate: true)
                .join(case_peaks, failOnMismatch: true, failOnDuplicate: true)
                .map { name, bams, bais, cps, ccs ->
                    def c = dcomps.find { x -> x.name == name }
                    tuple(name, c.ids.join(','), bams, bais,
                          c.ids.collect { sid -> sid in c.ctrlIds ? 'Control' : 'Condition' }.join(','),
                          c.blockCol ? c.ids.collect { sid -> c.blockOf[sid] }.join(',') : '',
                          cps, ccs, c.computed.join(','), c.fixedDisp ?: params.diff_fixed_dispersion)
                }
            DIFF_ACCESS(diff_in, txdb, params.genome_build, 'UCSC', params.diff_minq, params.diff_pe, params.diff_max_frag,
                        dn.diffDedup, params.diff_bin_size, params.diff_lfc, params.diff_fdr, params.diff_table,
                        params.diff_qlf_min_libs, params.diff_promoter_up, params.diff_promoter_down, params.diff_promoter_min_fc,
                        params.plot_font, params.plot_format, sc.diff_accessibility)
            diff_out = DIFF_ACCESS.out.peaks.map { n, f -> tuple(n, 'peaks', f) }
                .mix(DIFF_ACCESS.out.promoters.map { n, f -> tuple(n, 'promoters', f) })
                .mix(diff_out)
        }

        if (dn.runAnn) {
            // One task per comparison, on its peaks table (computed, or from diff_tables).
            DIFF_ANNOTATE(diff_out.filter { n, r, f -> r == 'peaks' }.map { n, r, f -> tuple(n, f) }, txdb, params.species,
                          params.genome_build, dn.bsgenome, dn.orgDb, params.annotation_tss_region, textOf(params.annotation_fdr),
                          params.plot_font, params.plot_format, sc.diff_annotate)
            // The DIFF_ANNOTATE tables are named NAME_<set>.csv: the set is the file name without the comparison prefix.
            annotation_out = DIFF_ANNOTATE.out.tables.flatMap { n, files ->
                [files].flatten().collect { f -> tuple(n, f.name.substring(n.length() + 1).replaceFirst(/\.csv$/, ''), f) }
            }
        }

        if (dn.runMot) {
            // Diff mode: one task per (comparison, motif region) with the matching table (differential or annotation).
            def motif_tables = diff_out.mix(annotation_out).filter { n, r, f -> r in dn.motifRegions && (n in dn.compNames) }
            def motif_in = motif_tables.map { n, r, f ->
                tuple(n, r, 'diff', f, '', motifLabelOf(r), params.motif_min_bin_peaks, no_file)
            }
            // Target mode: one task per motif_targets row, on the all_sites table of its comparison (the universe).
            if (dn.motifTargets) {
                def all_sites = annotation_out.filter { n, s, f -> s == 'all_sites' }.map { n, s, f -> tuple(n, f) }
                def targets_ch = Channel.fromList(dn.motifTargets.collect { t -> tuple(t.comparison, t) })
                def target_in = targets_ch.combine(all_sites, by: 0).map { comp, t, f ->
                    tuple(t.name, '', t.type, f, t.type == 'gene' ? t.target : '', t.label, t.minPeaks, t.bed ?: no_file)
                }
                motif_in = motif_in.mix(target_in)
            }
            MOTIF_ENRICHMENT(motif_in, params.species, params.genome_build, dn.bsgenome, textOf(params.motif_pwm_rds),
                             params.motif_min_score, params.motif_test, params.motif_max_kmer, params.motif_background,
                             params.motif_seed, params.motif_top, params.motif_padj_primary, params.motif_padj_relaxed,
                             params.motif_lfc, params.motif_gene_column, params.motif_match, params.plot_font, params.plot_format,
                             sc.motif_enrichment)

            // MOTIF_COMMON per region set, over the NON-pooled comparisons with at least two tables of that set. The tables and
            // summaries are grouped by region with groupTuple() and no size. A skipped analysis (NAME.skipped.txt only) emits
            // nothing on these tuples, so no group waits for a count that a skip never reaches; the groups are released when
            // MOTIF_ENRICHMENT has finished. The names are put in comparison order.
            def common_tables = MOTIF_ENRICHMENT.out.table
                .filter { n, r, f -> r in dn.motifRegions && (n in dn.compNames) && !(n in dn.pooledNames) }
                .map { n, r, f -> tuple(r, n, f) }
                .groupTuple()
                .filter { r, ns, fs -> ns.size() >= 2 }
            def common_summaries = MOTIF_ENRICHMENT.out.summary
                .filter { n, r, f -> r in dn.motifRegions && (n in dn.compNames) && !(n in dn.pooledNames) }
                .map { n, r, f -> tuple(r, n, f) }
                .groupTuple()
            def common_in = common_tables.join(common_summaries, failOnMismatch: true, failOnDuplicate: true)
                .map { r, ns, fs, ns2, ss ->
                    def order = (0..<ns.size()).toList().sort { i -> dn.compNames.indexOf(ns[i]) }
                    def order2 = (0..<ns2.size()).toList().sort { i -> dn.compNames.indexOf(ns2[i]) }
                    tuple(r, order.collect { i -> ns[i] }.join(','), order.collect { i -> fs[i] },
                          order2.isEmpty() ? [no_file] : order2.collect { i -> ss[i] })
                }
            MOTIF_COMMON(common_in, textOf(params.motif_common_label_strip), params.motif_common_padj, params.motif_common_min_lines,
                         params.plot_font, params.plot_format, sc.motif_common)
        }
    }

    // ---- F. MultiQC ------------------------------------------------------------------
    if (isTrue(params.run_multiqc)) {
        def multiqc_in = fastqc_zips.mix(trim_logs, align_logs, dedup_logs, raw_flagstat,
                                         raw_idxstats, final_flagstat, final_idxstats,
                                         macs_xls).collect()
        MULTIQC(multiqc_in, sc.multiqc)
    }
}
