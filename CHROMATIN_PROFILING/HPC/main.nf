#!/usr/bin/env nextflow
/*
 * CUTANDRUN - reproducible CUT&RUN / CUT&Tag pipeline
 * ======================================================
 * FASTQ (or pre-aligned BAM) -> trim -> align (target + spike-in genome) -> filter ->
 * dedup / spike split -> blacklist -> QC -> normalisation (ChIPseqSpikeInFree | spike-in |
 * TMM | CPM) -> bigWigs -> csaw peaks (default), MACS3 / SEACR (optional) -> QC report.
 * Optional: heatmaps, genome tracks, peak annotation, peak-set intersections.
 *
 * Run from your own project directory with:
 *   nextflow run /project/rephilllab/CUTANDRUN -profile pmacs -params-file config.yaml
 * or via the bundled wrapper:  run_cutandrun.sh
 *
 * Inputs (project root): config.yaml, samplesheet.csv, norm_runs.csv, peaksets.csv
 */

nextflow.enable.dsl = 2

include { ALIGN } from './modules/local/align'
include { BAM_FILTER } from './modules/local/bam_filter'
include { BAM_FINALIZE } from './modules/local/bam_finalize'
include { BAM_PREALIGNED } from './modules/local/bam_prealigned'
include { BIGWIG_LOG2RATIO } from './modules/local/bigwig_log2ratio'
include { BIN_COUNTS } from './modules/local/bin_counts'
include { BIGWIG } from './modules/local/bigwig'
include { BLACKLIST_FILTER } from './modules/local/blacklist_filter'
include { BUILD_MATRIX } from './modules/local/heatmap_build_matrix'
include { CHIPSEQSPIKEINFREE } from './modules/local/chipseqspikeinfree'
include { CLASSIFY_REGIONS } from './modules/local/classify_regions'
include { CLASSIFY_REGIONS_PLOT } from './modules/local/classify_regions_plot'
include { COMPUTEMATRIX } from './modules/local/computematrix'
include { COMPUTEMATRIX as COMPUTEMATRIX_CLUSTERED } from './modules/local/computematrix'
include { CORRELATION } from './modules/local/correlation'
include { CSAW_CUTOFF_SCAN } from './modules/local/csaw_cutoff_scan'
include { CSAW } from './modules/local/csaw'
include { DEDUP_SPLIT } from './modules/local/dedup_split'
include { FASTQC } from './modules/local/fastqc'
include { FINGERPRINT } from './modules/local/fingerprint'
include { FRIP } from './modules/local/frip'
include { HEATMAP_CLUSTER } from './modules/local/heatmap_cluster'
include { HEATMAP_LIMITS } from './modules/local/heatmap_limits'
include { INTERSECT } from './modules/local/intersect'
include { LIBCOMPLEXITY_FRAGS } from './modules/local/libcomplexity_frags'
include { MACS3 } from './modules/local/macs3'
include { MARKDUP } from './modules/local/markdup'
include { MULTIQC } from './modules/local/multiqc'
include { NORM_COMBINE } from './modules/local/norm_combine'
include { NORM_CPM } from './modules/local/norm_cpm'
include { NORM_QC } from './modules/local/norm_qc'
include { NORM_QC_PLOT } from './modules/local/norm_qc_plot'
include { PEAK_ANNOTATION } from './modules/local/peak_annotation'
include { PEAK_QC_REGIONS } from './modules/local/peak_qc_regions'
include { PICARD_ELC } from './modules/local/picard_elc'
include { PLOT_ENRICHED } from './modules/local/plot_enriched'
include { PLOTHEATMAP } from './modules/local/plotheatmap'
include { PRESEQ } from './modules/local/preseq'
include { PYGENOMETRACKS } from './modules/local/pygenometracks'
include { QC_ALIGNMENT_PLOT } from './modules/local/qc_alignment_plot'
include { QC_ALIGNMENT_TABLE } from './modules/local/qc_alignment_table'
include { QC_FRAGMENT_SIZES_PLOT } from './modules/local/qc_fragment_sizes_plot'
include { QC_FRAGMENT_SIZES_TABLE } from './modules/local/qc_fragment_sizes_table'
include { QC_FRIP_TABLE } from './modules/local/qc_frip_table'
include { QC_LIBCOMPLEXITY_PLOT } from './modules/local/qc_libcomplexity_plot'
include { QC_LIBCOMPLEXITY_TABLE } from './modules/local/qc_libcomplexity_table'
include { QC_PEAK_STATS } from './modules/local/qc_peak_stats'
include { QC_PRESEQ_PLOT } from './modules/local/qc_preseq_plot'
include { QC_TSS_TABLE } from './modules/local/qc_tss_table'
include { REPLICATE_CORRELATION } from './modules/local/replicate_correlation'
include { REPLICATE_CORRELATION_PLOT } from './modules/local/replicate_correlation_plot'
include { REPORT } from './modules/local/report'
include { SEACR_BEDGRAPH } from './modules/local/seacr_bedgraph'
include { SEACR_CALL } from './modules/local/seacr_call'
include { SEACR_CONSENSUS } from './modules/local/seacr_consensus'
include { SPIKEIN } from './modules/local/spikein'
include { TMM } from './modules/local/tmm'
include { TRACKS_INI } from './modules/local/tracks_ini'
include { TRACKS_LIMITS } from './modules/local/tracks_limits'
include { TRIM } from './modules/local/trim'
include { isBroadMark } from './modules/local/functions'
include { TSS_ENRICHMENT } from './modules/local/tss_enrichment'

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
    def rows = System.getProperties().get('cutandrun.host_paths')
    if (rows == null) {
        rows = []
        System.getProperties().put('cutandrun.host_paths', rows)
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

// call_peaks (default true): false switches off every peak caller and everything that only consumes called peaks.
def callPeaks() { return params.call_peaks == null || isTrue(params.call_peaks) }
def peakOn(String k) { return callPeaks() && isTrue(params[k]) }

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

// The six threshold columns of a csaw peak set (peaksets.csv order).
def csawThresholdCols() {
    return ['small_filt', 'large_filt', 'merge_bp', 'max_width', 'fdr', 'lfc']
}

// csaw_defaults entry of one antibody (exact key, else case-insensitive); empty map when there is none.
def csawDefaultsFor(String antibody) {
    def defs = params.csaw_defaults ?: [:]
    if (!antibody || !(defs instanceof Map)) return [:]
    def key = defs.containsKey(antibody) ? antibody : defs.keySet().find { k -> k.toString().equalsIgnoreCase(antibody) }
    if (key == null || !(defs[key] instanceof Map)) return [:]
    return defs[key]
}

// Every csaw_defaults entry must be a map whose keys are threshold columns; returns the messages.
def csawDefaultsProblems() {
    def defs = params.csaw_defaults ?: [:]
    if (!(defs instanceof Map)) return []   // reported by prepareRun's parameter checks
    def out = []
    defs.each { ab, entry ->
        if (!(entry instanceof Map)) {
            out << "csaw_defaults.${ab} must be a map of ${csawThresholdCols().join(', ')}"
        } else {
            def bad = entry.keySet().findAll { k -> !(k.toString() in csawThresholdCols()) }
            if (bad) out << "csaw_defaults.${ab}: unknown key(s) ${bad.join(', ')}; allowed keys: ${csawThresholdCols().join(', ')}"
        }
    }
    return out
}

// Writes the resolved csaw thresholds (value and source of each) to metadata/peaksets.resolved.tsv.
// A record only: a failure is a warning, not a stop.
def writePeaksetsResolved(List prows) {
    def cols = csawThresholdCols()
    def header = ['peakset', 'target_group', 'control_group', 'antibody'] + cols + cols.collect { c -> "${c}_source".toString() }
    try {
        def out = file("${params.metadata_dir}/peaksets.resolved.tsv")
        out.parent.mkdirs()
        def lines = [header.join('\t')]
        prows.each { p ->
            lines << ([p.peakset, p.target_group, p.control_group, p.antibody] + cols.collect { c -> p[c] } + cols.collect { c -> p.sources[c] }).join('\t')
        }
        out.text = lines.join('\n') + '\n'
    } catch (Exception e) {
        log.warn "could not write ${params.metadata_dir}/peaksets.resolved.tsv: ${e.message}"
    }
}

// Peak-calling mode of every target group, from its antibody (broad_marks): MACS3 broad / narrow, SEACR relaxed / stringent.
// macs3_broad = true | false and seacr_mode = relaxed | stringent force every group; auto follows isBroadMark.
// Returns [rows: [group, antibody, caller, mode, rule], macsBroad: group -> boolean, seacrMode: sample_id -> mode, errs: [...]].
def peakCallModes(Map groups, Set controlIds, Map controlOf, Map meta, boolean runMacs3, boolean runSeacr) {
    def macsForce = textOf(params.macs3_broad).toLowerCase()
    def seacrForce = textOf(params.seacr_mode).toLowerCase()
    def rows = []
    def errs = []
    def macsBroad = [:]
    def seacrMode = [:]
    groups.each { g, ids ->
        def targets = ids.findAll { id -> !controlIds.contains(id) }
        def called = targets.findAll { id -> controlOf.containsKey(id) }
        [[caller: 'MACS3', on: runMacs3, members: targets, forced: macsForce in ['true', 'false']],
         [caller: 'SEACR', on: runSeacr, members: called, forced: seacrForce in ['relaxed', 'stringent']]].each { c ->
            if (c.on && c.members) {
                def abs = c.members.collect { id -> meta[id].antibody.toString() }.unique()
                def flags = abs.collect { ab -> isBroadMark(ab, params.broad_marks) }.unique()
                def broad
                if (c.forced) {
                    broad = c.caller == 'MACS3' ? macsForce == 'true' : seacrForce == 'relaxed'
                } else {
                    if (flags.size() > 1) errs << "group ${g}: antibodies ${abs.join(', ')} differ in broad / narrow status (broad_marks: ${params.broad_marks}); ${c.caller} calls one mode per group. Put them in separate groups, or force the mode with ${c.caller == 'MACS3' ? 'macs3_broad' : 'seacr_mode'}."
                    broad = flags[0]
                }
                def mode = c.caller == 'MACS3' ? (broad ? 'broad' : 'narrow') : (broad ? 'relaxed' : 'stringent')
                rows << [group: g, antibody: abs.join('_'), caller: c.caller, mode: mode, rule: c.forced ? 'forced' : 'auto:broad_marks']
                if (c.caller == 'MACS3') macsBroad[g] = broad
                else called.each { id -> seacrMode[id] = mode }
            }
        }
    }
    return [rows: rows, macsBroad: macsBroad, seacrMode: seacrMode, errs: errs]
}

// Writes the peak-calling mode of every group and caller to metadata/peak_call_modes.tsv. A record only: a failure is a warning.
def writePeakCallModes(List rows) {
    try {
        def out = file("${params.metadata_dir}/peak_call_modes.tsv")
        out.parent.mkdirs()
        def lines = [['group', 'antibody', 'caller', 'mode', 'rule'].join('\t')]
        rows.each { r -> lines << [r.group, r.antibody, r.caller, r.mode, r.rule].join('\t') }
        out.text = lines.join('\n') + '\n'
    } catch (Exception e) {
        log.warn "could not write ${params.metadata_dir}/peak_call_modes.tsv: ${e.message}"
    }
}

// Text listing the csaw peak sets, the classification regions and the intersect names a region reference may name.
def regionDesc(List csawNames, List classNames, List intersectNames) {
    return "csaw peak sets: ${csawNames ? csawNames.join(', ') : (callPeaks() ? 'none' : 'none (call_peaks is false: no peak sets are called)')}; classification regions: ${classNames ? classNames.join(', ') : 'none'}; intersections: ${intersectNames ? intersectNames.join(', ') : 'none'}"
}

// Comma list of a config or CSV value, trimmed; blank entries are kept (null or blank = []).
def csvListKeep(v) {
    if (v == null || v.toString().trim() == '') return []
    return v.toString().split(',', -1).collect { s -> s.trim() }
}

// A region / BED reference: a csaw peak set called in this run ([kind: 'peakset']), a classification
// region ([kind: 'class']), an intersect name ([kind: 'intersect']) or an existing BED file ([kind: 'path']).
// Returns null when it is none of these.
def classifyBed(String ref, List csawNames, List classNames, List intersectNames) {
    if (ref in csawNames) return [kind: 'peakset', value: ref]
    if (ref in classNames) return [kind: 'class', value: ref]
    if (ref in intersectNames) return [kind: 'intersect', value: ref]
    def hp = hostPath(ref)   // a Windows drive path is a '/' path inside the container
    def bedLike = hp.contains('/') || hp.toLowerCase().endsWith('.bed') || hp.toLowerCase().endsWith('.bed.gz')
    if (bedLike && file(hostPathNote('region or BED reference', ref)).exists()) return [kind: 'path', value: ref]
    return null
}

// Default label of a BED reference: the peak set / classification region name, or the file name without its extension.
def bedLabel(String ref, Map bm) {
    if (bm == null || bm.kind in ['peakset', 'class', 'intersect']) return ref
    def nm = file(hostPath(ref)).name.replaceFirst('(?i)\\.gz$', '')
    return nm.contains('.') ? nm.substring(0, nm.lastIndexOf('.')) : nm
}

// labels / colors: empty (default) or exactly one per sample.
def countProblems(String lbl, String col, List vals, int n) {
    if (vals && vals.size() != n) return ["${lbl}: ${col} has ${vals.size()} entries but samples has ${n}"]
    return []
}

// Removes the entries of one list column of a CSV row for which gone(entry) is true, and the same positions in the
// columns aligned with it (a column whose entry count differs, e.g. one 'auto' for all samples, is left alone).
// Returns [row: the row (a copy when something was removed), dropped: the removed entries].
def dropFromRow(Map row, String keyCol, gone, List aligned) {
    def keys = csvListKeep(row[keyCol])
    def idx = (0..<keys.size()).findAll { i -> !gone.call(keys[i]) }
    if (idx.size() == keys.size()) return [row: row, dropped: []]
    def out = new LinkedHashMap(row)
    out[keyCol] = idx.collect { i -> keys[i] }.join(',')
    aligned.each { col ->
        def parts = csvListKeep(row[col])
        if (parts.size() == keys.size()) out[col] = idx.collect { i -> parts[i] }.join(',')
    }
    return [row: out, dropped: keys.findAll { k -> gone.call(k) }.unique()]
}

// A heatmaps.csv / tracks.csv row against the samplesheet: samples that are not active (absent from the samplesheet,
// excluded or in skip_groups) are dropped from the row with their positions in the per-sample columns, regions / beds
// that name a skipped peak set, classification or intersect are dropped too. Prints one WARN per change. Returns the
// pruned row, or null (with a WARN) when no sample is left (or, with needRegions, no region).
def pruneRow(String csv, String sheet, String what, Map row, List activeIds, List skipped, List sampleCols, String regionCol, List regionCols, boolean needRegions) {
    def name = "${what} ${row[what]}"
    def sr = dropFromRow(row, 'samples', { x -> x && !(x in activeIds) }, sampleCols)
    def cur = sr.row
    if (sr.dropped) {
        log.warn "${csv}: ${name}: dropped sample(s) not in ${sheet} or not active (excluded / skip_groups): ${sr.dropped.join(', ')}"
        if (!csvList(cur.samples)) {
            log.warn "${csv}: ${name} skipped: none of its samples is in ${sheet}"
            return null
        }
        if (cur.containsKey('sort_using') && cur.sort_using in sr.dropped) {
            log.warn "${csv}: ${name}: sort_using '${cur.sort_using}' was dropped; default sort used"
            cur = new LinkedHashMap(cur); cur.sort_using = ''
        }
        cur = dropFromRow(cur, 'cluster_samples', { x -> x in sr.dropped }, ['cluster_labels']).row
    }
    def rr = dropFromRow(cur, regionCol, { x -> x in skipped }, regionCols)
    if (rr.dropped) {
        log.warn "${csv}: ${name}: dropped ${regionCol} entry(ies) of skipped rows (no active samples): ${rr.dropped.join(', ')}"
        if (needRegions && !csvList(rr.row[regionCol])) {
            log.warn "${csv}: ${name} skipped: none of its ${regionCol} is left"
            return null
        }
    }
    return rr.row
}

// samples of a heatmap or track row: non-empty, active samples only, each listed once.
// csvPath = the CSV the row comes from, allIds = every sample_id of the samplesheet; a missing id is reported
// with the CSV, the row and the fix, all ids of the row in one message.
def sampleProblems(String lbl, List samples, List activeIds, String sheet, String csvPath = '', List allIds = []) {
    def out = []
    if (!samples) out << "${lbl}: samples is empty"
    def unknown = samples.findAll { s -> !(s in allIds) && !(s in activeIds) }.unique()
    def inactive = samples.findAll { s -> (s in allIds) && !(s in activeIds) }.unique()
    if (unknown) out << "${csvPath ?: sheet}: ${lbl}: sample(s) not in the samplesheet ${sheet}: ${unknown.join(', ')}. Remove them from ${csvPath ?: 'this CSV'} or add them to the samplesheet."
    if (inactive) out << "${csvPath ?: sheet}: ${lbl}: sample(s) in ${sheet} but not active (excluded or in skip_groups): ${inactive.join(', ')}. Remove them from ${csvPath ?: 'this CSV'} or make them active in the samplesheet."
    def twice = samples.findAll { s -> samples.count(s) > 1 }.unique()
    if (twice) out << "${lbl}: samples lists a sample twice: ${twice.join(', ')}"
    return out
}

// z_min, z_max, y_max, max_value: none, one, or one number per sample.
def numListProblems(String lbl, String col, List vals, int n) {
    if (!(vals.size() in [0, 1, n])) return ["${lbl}: ${col} has ${vals.size()} entries; give none, one, or one per sample (${n})"]
    return vals.findAll { v -> !isNum(v) }.collect { v -> "${lbl}: ${col} entry '${v}' is not a number" }
}

// z_max / y_max style lists: none, one, or one per sample; each a number or (allowAuto) 'auto'.
def limitListProblems(String lbl, String col, List vals, int n, boolean allowAuto) {
    if (!(vals.size() in [0, 1, n])) return ["${lbl}: ${col} has ${vals.size()} entries; give none, one, or one per sample (${n})"]
    def form = allowAuto ? 'a number or auto' : 'a number'
    return vals.findAll { v -> !(isNum(v) || (allowAuto && v == 'auto')) }.collect { v -> "${lbl}: ${col} entry '${v}' is not ${form}" }
}

// One max_value token per sample of a tracks.csv row: number | auto | group:NAME. A blank entry (or no entry) takes dflt;
// 'shared' becomes the group of the sample's antibody (group:shared_<antibody>, other characters replaced by _).
def trackMaxTokens(List entries, int n, String dflt, List antibodies) {
    return (0..<n).collect { i ->
        def e = entries.isEmpty() ? '' : (entries.size() == 1 ? entries[0] : entries[i])
        if (!e) e = dflt
        if (e == 'shared') return "group:shared_${(antibodies[i] ?: "unknown").toString().replaceAll('[^A-Za-z0-9_.-]', '_')}".toString()
        return e.toString()
    }
}

// ---- Peak QC (run_peak_qc): everything known before the run starts --------------------------------------------
def peakQcMark(Map cfg, List ids) {
    return ids.collect { id -> cfg.meta[id].antibody.toString() }.unique().join('_')
}

// A classification id or group name without its mark token: 'EZH2i500_K27me3' -> 'EZH2i500'.
def peakQcStrip(String s, String mark) {
    return s.split('_').findAll { t -> !t.equalsIgnoreCase(mark) }.join('_')
}

// Peak sets of this run (csaw: <peakset>; MACS3: macs3_<group>; SEACR: seacr_<group>) with their mark and samples, the
// classifications whose base set is a csaw peak set, the bigWig samples and the samples/colours/y-max tokens of the track set peak_qc.
def peakQcPlan(Map cfg) {
    def sets = [:]
    if (peakOn('run_csaw')) {
        cfg.peaksets.each { ps ->
            sets[ps.peakset.toString()] = [mark: peakQcMark(cfg, ps.target_ids), targets: ps.target_ids,
                controls: (ps.target_ids.collect { id -> cfg.controlOf[id] } + ps.control_ids).findAll { c -> c }.unique()]
        }
    }
    def tgroups = cfg.groups.collect { g, ids -> [group: g, targets: ids.findAll { id -> !cfg.isControl.contains(id) }] }.findAll { gr -> gr.targets }
    if (cfg.runMacs3) {
        tgroups.each { gr ->
            sets["macs3_${gr.group}".toString()] = [mark: peakQcMark(cfg, gr.targets), targets: gr.targets,
                controls: gr.targets.collect { id -> cfg.controlOf[id] }.findAll { c -> c }.unique()]
        }
    }
    if (cfg.runSeacr && cfg.controlOf) {
        tgroups.each { gr ->
            def called = gr.targets.findAll { id -> cfg.controlOf.containsKey(id) }
            if (called) {
                sets["seacr_${gr.group}".toString()] = [mark: peakQcMark(cfg, called), targets: called,
                    controls: called.collect { id -> cfg.controlOf[id] }.findAll { c -> c }.unique()]
            }
        }
    }
    def marks = sets.values().collect { s -> s.mark }.unique()

    def classes = [:]
    def clsList = []
    if (isTrue(params.run_classify) && peakOn('run_csaw')) {
        cfg.classifications.each { c ->
            if (c.regions.kind == 'peakset' && sets.containsKey(c.regions.value.toString())) {
                def mark = sets[c.regions.value.toString()].mark
                def contrast = peakQcStrip(c.classification.toString(), mark)
                if (!contrast || contrast == c.classification.toString()) {
                    contrast = "${peakQcStrip(c.ref_group.toString(), mark)}_vs_${peakQcStrip(c.treat_group.toString(), mark)}".toString()
                }
                classes[c.classification.toString()] = [mark: mark, base: c.regions.value.toString(), contrast: contrast, ref_group: c.ref_group.toString()]
                clsList << c
            }
        }
    }

    def bwIds = (sets.values().collectMany { s -> s.targets } + sets.values().collectMany { s -> s.controls }).unique()
    // track samples: targets of the peak-set groups and of the classification groups, then their IgG controls
    def trackTargets = (sets.values().collectMany { s -> s.targets } + clsList.collectMany { c -> c.ref_ids + c.treat_ids }).unique()
    def trackControls = (sets.values().collectMany { s -> s.controls } + trackTargets.collect { id -> cfg.controlOf[id] }.findAll { c -> c }).unique()
    def trackIds = trackTargets.sort(false) { id -> "${cfg.meta[id].antibody}\t${cfg.meta[id].group}\t${id}".toString() } +
        trackControls.sort(false) { id -> "${cfg.meta[id].group}\t${id}".toString() }
    def palette = ['#8c6bb1', '#8c96c6', '#41ab5d', '#e08214', '#2a78d6', '#e34948', '#17becf', '#bcbd22']
    def trackAbs = trackIds.collect { id -> cfg.meta[id].antibody.toString() }
    def distinct = trackAbs.findAll { ab -> !ab.equalsIgnoreCase('IgG') }.unique()
    def colors = trackAbs.collect { ab -> ab.equalsIgnoreCase('IgG') ? '#808080' : palette[distinct.indexOf(ab) % palette.size()] }
    def maxTok = trackMaxTokens([], trackIds.size(), textOf(params.tracks_max_value_default) ?: 'shared', trackAbs)
    return [sets: sets, marks: marks, classes: classes, bw_ids: bwIds, track_ids: trackIds, track_colors: colors, track_max: maxTok,
            mark_of: cfg.meta.collectEntries { id, r -> [(id): r.antibody.toString()] },
            group_of: cfg.meta.collectEntries { id, r -> [(id): r.group.toString()] }]
}

// Colours must be #RRGGBB or a colour name (colorRe is the pattern).
def colourProblems(String lbl, List colors, String colorRe) {
    return colors.findAll { c -> !(c ==~ colorRe) }.collect { c -> "${lbl}: colour '${c}' is not a #RRGGBB hex or a colour name" }
}

def hasUnsafeChar(String s) {
    return s != null && (s.contains('"') || s.contains("'") || s.contains('`') || s.contains('$'))
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

def prepareRun() {
    if (params.help) {
        log.info """
        CUTANDRUN pipeline (CUT&RUN / CUT&Tag, paired-end)
        ---------------------------------------------------
        Usage:
          nextflow run ${workflow.projectDir} -profile pmacs -params-file config.yaml

        Required params (set in config.yaml):
          samplesheet     CSV: sample_id,fastq_1,fastq_2,bam,antibody,condition,replicate,control_id,group,norm_run,exclude
                          (fastq_1+fastq_2 OR bam per row; control = antibody IgG)
          genome_index    bowtie2 index prefix of the combined target + spike-in genome (FASTQ rows)
          chrom_sizes     chrom.sizes of the combined genome
          blacklist       blacklist BED (FASTQ rows)
          gtf             gene annotation GTF; TSS enrichment anchors on its transcript start sites (also needed by run_tracks,
                          run_annotation and run_peak_qc)
          tss_bed         TSS BED, fallback for TSS enrichment; required only when gtf is not set
          norm_runs       CSV: norm_run,method,bin_size,cutoff,max_turns (unless normalisation: cpm); method = the PRIMARY
                          method (chipseqspikeinfree, spikein, tmm or cpm); the other norm_methods are computed for comparison
                          optional columns: min_first_turn (chipseqspikeinfree; blank = 0.5, the package default),
                          ref_condition (reference condition of the normalisation-QC MA plots)
          peaksets        CSV: peakset,target_group,control_group,small_filt,large_filt,merge_bp,max_width,fdr,lfc
                          (when run_csaw or run_csaw_cutoff_scan; a blank threshold cell takes csaw_defaults[antibody],
                          then csaw_<col>; see nextflow.config)

        Optional CSV inputs (each needs its run_* switch and its path set in config.yaml):
          classifications CSV: classification,regions,ref_group,treat_group,norm_run,lfc_retained,lfc_lost,fdr
                          (run_classify; regions = a csaw peak set name or a BED path; blank norm_run = the one norm_run of
                          the samples; blank thresholds take classify_lfc_retained, classify_lfc_lost, classify_fdr)
          heatmaps       CSV: heatmap,regions,samples,labels,window_bp,bin_bp,sort_using,cluster_k,z_min,z_max,y_max,colors,engine
                          optional columns: cluster_samples,cluster_labels,legend_title (run_heatmaps)
          tracks          CSV: track_set,locus_or_gene,flank,samples,labels,colors,max_value,beds,scalebar
                          optional columns: log2ratio, bed_colors, bed_labels, ini_args (run_tracks; needs gtf; ini_args = extra make_tracks_ini.py options for that set)
          intersects      CSV: name,bed_a,bed_b,label_a,label_b,min_overlap,reciprocal[,output] (run_intersect)

        Optional: project_name, outdir (./results), experiment (cutandrun|cutandtag), skip_groups,
                  call_peaks (default true; false = no peak calling, FRiP, peak QC or peak annotation),
                  run_macs3, run_seacr, run_log2ratio, run_csaw, run_classify (needs classifications), run_report,
                  run_heatmaps, run_tracks, run_annotation (needs gtf and a peak source), run_intersect,
                  run_peak_qc (default true: automatic peak QC tracks; peak_qc_* options),
                  and the other keys in examples/config.yaml.
                  csaw knobs: csaw_small_filt, csaw_large_filt, csaw_merge_bp, csaw_max_width, csaw_fdr, csaw_lfc (fallbacks
                  for blank peaksets.csv cells), csaw_defaults (per-antibody fallbacks), csaw_large_window, csaw_large_spacing,
                  csaw_bg_max_enrich, csaw_count_filter, csaw_max_frag, csaw_lfc_from (best|combined), csaw_disjoint_regions. Replace a step script with  <script name>: path  in config.yaml.
        """
        System.exit(0)
    }

    def errs = []
    def idRe = '[A-Za-z0-9][A-Za-z0-9_.-]*'

    // ---- Top-level params ---------------------------------------------------
    if (!isTrue(params.paired_end)) fail("paired_end must be true: this pipeline version supports paired-end data only.")
    if (!params.samplesheet) fail("missing required param: samplesheet (set in your config.yaml).")
    if (params.publish_mode == 'auto') {
        fail("publish_mode 'auto' is resolved by run_cutandrun (it passes --publish_mode link or copy); when running Nextflow by hand pass --publish_mode copy")
    }

    // Enumerated params.
    def allowed = [bowtie2_mode: ['local', 'end_to_end'], experiment: ['cutandrun', 'cutandtag'],
                   bigwig_units: ['rpkm', 'cpm', 'none'], csaw_norm_bins: ['background', 'all'],
                   normalisation: ['chipseqspikeinfree', 'spikein', 'tmm', 'cpm'], seacr_mode: ['auto', 'relaxed', 'stringent'], macs3_broad: ['auto', 'true', 'false'],
                   csaw_lfc_from: ['best', 'combined'], tracks_limit_stat: ['mean_bins', 'max_bin']]
    allowed.each { k, vals ->
        def v = params[k] == null ? null : params[k].toString()
        if (!(v in vals)) errs << "${k}: '${v}' is not valid (one of: ${vals.join(', ')})"
    }
    // Numeric csaw knobs, and the optional max.frag override (null = bowtie2_maxins, else 500).
    ['csaw_large_window', 'csaw_large_spacing', 'csaw_bg_max_enrich', 'csaw_count_filter'].each { k ->
        if (!isNum(params[k])) errs << "${k}: '${params[k]}' must be a number"
    }
    if (textOf(params.csaw_max_frag)) {
        if (!isNum(params.csaw_max_frag) || params.csaw_max_frag.toString().toDouble() <= 0) errs << "csaw_max_frag: '${params.csaw_max_frag}' must be a positive number (or null for bowtie2_maxins)"
    }
    if (!(textOf(params.csaw_disjoint_regions).toLowerCase() in ['', 'true', 'false', '1', '0', 'yes', 'no'])) {
        errs << "csaw_disjoint_regions: '${params.csaw_disjoint_regions}' must be true or false"
    }
    if (params.csaw_defaults != null && !(params.csaw_defaults instanceof Map)) errs << "csaw_defaults must be a map of antibody -> {small_filt, large_filt, merge_bp, max_width, fdr, lfc}"
    ['classify_lfc_retained', 'classify_lfc_lost', 'classify_fdr'].each { k ->
        if (!isNum(params[k])) errs << "${k}: '${params[k]}' must be a number"
    }
    // Normalisation methods computed for every norm_run besides its primary, and the numeric normalisation-QC / replicate knobs.
    def normMethods = csvList(params.norm_methods).collect { m -> m.toLowerCase() }.unique()
    def badMethods = normMethods.findAll { m -> !(m in ['chipseqspikeinfree', 'tmm', 'cpm']) }
    if (badMethods) errs << "norm_methods: '${badMethods.join(', ')}' is not valid (comma list of chipseqspikeinfree, tmm, cpm; '' = primary method only; spike-in is switched on by norm_spikein)"
    ['csif_bin_size', 'csif_cutoff', 'csif_max_turns', 'norm_qc_min_count', 'norm_qc_prior', 'norm_qc_panels_per_page', 'repcor_prior', 'repcor_signal_frac'].each { k ->
        if (!isNum(params[k])) errs << "${k}: '${params[k]}' must be a number"
    }
    if (textOf(params.repcor_bin_size) && !isNum(params.repcor_bin_size)) errs << "repcor_bin_size: '${params.repcor_bin_size}' must be a number (or null for tmm_bin_size)"
    failList("invalid config.yaml value(s)", errs)
    def adapterSet = params.experiment.toString() == 'cutandtag' ? 'NexteraPE-PE' : 'TruSeq3-PE-2'

    // ---- Samplesheet ----------------------------------------------------------
    def sheetCols = ['sample_id', 'fastq_1', 'fastq_2', 'bam', 'antibody', 'condition', 'replicate', 'control_id', 'group', 'norm_run', 'exclude']
    def allRows = readTable('samplesheet', params.samplesheet, sheetCols)
    if (!allRows) fail("${params.samplesheet} has no data rows.")

    // sample_id: present, safe in file names, unique.
    def allIds = allRows.collect { row -> row.sample_id }
    allRows.eachWithIndex { row, i ->
        if (!row.sample_id) errs << "row ${i + 1} has an empty sample_id"
        else if (!(row.sample_id ==~ idRe)) errs << "sample_id '${row.sample_id}' must match ^${idRe}\$ (letters, digits, _ . -; it is used verbatim in file names)"
    }
    def dupIds = allIds.findAll { s -> s && allIds.count(s) > 1 }.unique()
    if (dupIds) errs << "duplicate sample_id values: ${dupIds.join(', ')}"
    // exclude flag.
    allRows.each { row ->
        def e = row.exclude.toLowerCase()
        if (e in ['', 'false', '0', 'no']) row.excluded = false
        else if (e in ['true', '1', 'yes']) row.excluded = true
        else errs << "${row.sample_id}: exclude '${row.exclude}' must be empty/false/0/no or true/1/yes"
    }
    failList("${params.samplesheet}: invalid sample_id / exclude", errs)

    def allById = allRows.collectEntries { row -> [(row.sample_id): row] }
    def rows = allRows.findAll { row -> !row.excluded }   // everything below ignores excluded rows
    def excludedIds = allRows.findAll { row -> row.excluded }.collect { row -> row.sample_id }
    if (!rows) fail("${params.samplesheet}: every row is excluded.")

    // antibody, group.
    rows.each { row ->
        if (!row.antibody) errs << "${row.sample_id}: antibody is empty"
        if (!row.group) {
            if (row.condition && row.antibody) row.group = "${row.condition}_${row.antibody}".toString()
            else errs << "${row.sample_id}: group is empty and cannot be derived (needs condition and antibody)"
        }
        if (row.group && !(row.group ==~ idRe)) errs << "${row.sample_id}: group '${row.group}' must match ^${idRe}\$ (it is used in file names)"
    }
    failList("${params.samplesheet}: invalid antibody / group", errs)

    // Input files: exactly one of (fastq_1 AND fastq_2) or bam.
    rows.each { row ->
        def hasFq = row.fastq_1 && row.fastq_2
        if (row.bam && (row.fastq_1 || row.fastq_2)) {
            errs << "${row.sample_id}: give either fastq_1+fastq_2 or bam, not both"
        } else if (row.bam) {
            if (!row.bam.endsWith('.bam')) errs << "${row.sample_id}: bam must end in .bam (${row.bam})"
            if (!file(hostPathNote("samplesheet ${row.sample_id} bam", row.bam)).exists()) errs << "${row.sample_id}: bam not found: ${row.bam}"
        } else if (hasFq) {
            if (row.fastq_1 == row.fastq_2) errs << "${row.sample_id}: fastq_1 and fastq_2 are the same file (${row.fastq_1})"
            if (!file(hostPathNote("samplesheet ${row.sample_id} fastq_1", row.fastq_1)).exists()) errs << "${row.sample_id}: fastq_1 not found: ${row.fastq_1}"
            if (!file(hostPathNote("samplesheet ${row.sample_id} fastq_2", row.fastq_2)).exists()) errs << "${row.sample_id}: fastq_2 not found: ${row.fastq_2}"
        } else {
            errs << "${row.sample_id}: needs fastq_1 AND fastq_2 (paired-end), or bam"
        }
    }
    failList("${params.samplesheet}: input file problems", errs)

    // Skipped groups: kept (and validated above) but out of every core step.
    def skipGroups = (params.skip_groups ?: '').toString().split(',').collect { g -> g.trim() }.findAll { g -> g }
    def rowGroups = rows.collect { row -> row.group }.unique()
    def unusedSkips = skipGroups.findAll { g -> !(g in rowGroups) }
    if (unusedSkips) log.warn "skip_groups names no group in ${params.samplesheet}: ${unusedSkips.join(', ')}"
    def skippedIds = rows.findAll { row -> row.group in skipGroups }.collect { row -> row.sample_id }
    def activeRows = rows.findAll { row -> !(row.sample_id in skippedIds) }
    def activeIds = activeRows.collect { row -> row.sample_id }
    if (!activeRows) fail("every non-excluded row is in skip_groups (${skipGroups.join(', ')}): nothing to run.")

    // Controls (antibody IgG) and the target -> control links.
    def controlIds = rows.findAll { row -> row.antibody.equalsIgnoreCase('IgG') }.collect { row -> row.sample_id } as Set
    // call_peaks: false = no peak caller runs (MACS3 and SEACR are forced off; an explicit true is an error).
    if (!callPeaks()) {
        def explicitCallers = ['run_macs3', 'run_seacr'].findAll { k -> isTrue(params[k]) }
        if (explicitCallers) fail("call_peaks is false but ${explicitCallers.join(' and ')} ${explicitCallers.size() > 1 ? 'are' : 'is'} true: remove it or set call_peaks: true")
        log.info "call_peaks is false: no peaks are called (csaw, csaw cutoff scan, MACS3, SEACR, FRiP, peak QC and peak annotation are off; peaksets is ignored)"
    }
    // run_macs3 / run_seacr: not set (null) = on when every target has a control, else skipped with a warning; true / false = explicit.
    def macs3Auto = textOf(params.run_macs3) == ''
    def seacrAuto = textOf(params.run_seacr) == ''
    def runMacs3 = !callPeaks() ? false : (macs3Auto ? true : isTrue(params.run_macs3))
    def runSeacr = !callPeaks() ? false : (seacrAuto ? true : isTrue(params.run_seacr))
    def needsControl = (callPeaks() && ((!macs3Auto && runMacs3) || (!seacrAuto && runSeacr))) || isTrue(params.run_log2ratio)
    def noControl = []
    def controlOf = [:]
    rows.each { row ->
        def c = row.control_id
        def skippedTarget = row.sample_id in skippedIds
        if (row.sample_id in controlIds) {
            if (c) errs << "${row.sample_id}: a control (IgG) must have an empty control_id (has '${c}')"
        } else if (!c) {
            if (!skippedTarget) noControl << row.sample_id
        } else if (!allById.containsKey(c)) {
            errs << "${row.sample_id}: control_id '${c}' is not a sample_id in ${params.samplesheet}"
        } else if (allById[c].excluded) {
            errs << "${row.sample_id}: control_id ${c} is excluded"
        } else if (!(c in controlIds)) {
            errs << "${row.sample_id}: control_id '${c}' is not an IgG sample (antibody ${allById[c].antibody})"
        } else if ((c in skippedIds) && !skippedTarget) {
            errs << "${row.sample_id}: control_id '${c}' is in a skipped group but this sample is not"
        } else if (!skippedTarget) {
            controlOf[row.sample_id] = c
        }
    }
    failList("${params.samplesheet}: control_id problems", errs)
    if (noControl) {
        def msg = "no control_id for: ${noControl.join(', ')}"
        if (needsControl) fail("${msg}. MACS3 / SEACR / log2ratio need an IgG control per target (set control_id, or turn off run_macs3, run_seacr and run_log2ratio).")
        else {
            def skipped = []
            if (macs3Auto) { runMacs3 = false; if (callPeaks()) skipped << 'MACS3' }
            if (seacrAuto) { runSeacr = false; if (callPeaks()) skipped << 'SEACR' }
            if (skipped) log.warn "${msg}: ${skipped.join(' and ')} (on by default) are skipped; set control_id for every target to call them, or set run_macs3 / run_seacr to false to silence this. csaw is not affected."
            else if (callPeaks()) log.warn "${msg} (fine for csaw-only runs)"
        }
    }

    // Per-library tables.
    def fqById = [:]
    def bamById = [:]
    activeRows.each { row ->
        if (row.bam) bamById[row.sample_id] = file(hostPath(row.bam))
        else fqById[row.sample_id] = [file(hostPath(row.fastq_1)), file(hostPath(row.fastq_2))]
    }
    def fastqRows = fqById.collect { id, fq -> [id, fq[0], fq[1]] }
    def bamRows = bamById.collect { id, b -> [id, b] }
    def meta = rows.collectEntries { row -> [(row.sample_id): row] }
    def groups = [:]    // group -> active sample_ids, samplesheet order
    def groupsOf = [:]  // sample_id -> group name (all non-excluded rows)
    rows.each { row -> groupsOf[row.sample_id] = row.group }
    activeRows.each { row ->
        if (!groups.containsKey(row.group)) groups[row.group] = []
        groups[row.group] << row.sample_id
    }

    // ---- Peak-calling modes (MACS3 broad / narrow, SEACR relaxed / stringent) per group ----------
    def callModes = peakCallModes(groups, controlIds, controlOf, meta, runMacs3, runSeacr)
    failList("peak calling mode (broad_marks, macs3_broad, seacr_mode): problems", callModes.errs)
    callModes.rows.each { r -> log.info "peak calling: group ${r.group} (${r.antibody}) ${r.caller} ${r.mode} [${r.rule}]" }
    if (callModes.rows) writePeakCallModes(callModes.rows)

    // ---- Normalisation runs ---------------------------------------------------
    def normRuns = []
    def normRunOf = [:]   // sample_id -> norm_run name (active rows)
    if (params.normalisation.toString() != 'cpm') {
        if (!params.norm_runs) fail("norm_runs is required unless normalisation is cpm (set norm_runs in your config.yaml).")
        def nrows = readTable('norm_runs', params.norm_runs, ['norm_run', 'method', 'bin_size', 'cutoff', 'max_turns'])
        def names = nrows.collect { n -> n.norm_run }
        nrows.eachWithIndex { n, i ->
            if (!n.norm_run) {
                errs << "norm_runs row ${i + 1} has an empty norm_run"
            } else {
                if (!(n.norm_run ==~ idRe)) errs << "norm_run '${n.norm_run}' must match ^${idRe}\$ (it is used in file names)"
                def m = n.method.toLowerCase()
                if (!(m in ['chipseqspikeinfree', 'spikein', 'tmm', 'cpm'])) {
                    errs << "norm_run ${n.norm_run}: method '${n.method}' must be chipseqspikeinfree, spikein, tmm or cpm"
                } else if (m == 'chipseqspikeinfree') {
                    ['bin_size', 'cutoff', 'max_turns'].each { col ->
                        if (!isNum(n[col])) errs << "norm_run ${n.norm_run}: ${col} '${n[col]}' must be numeric for chipseqspikeinfree"
                    }
                    if (textOf(n.min_first_turn) && !isNum(n.min_first_turn)) errs << "norm_run ${n.norm_run}: min_first_turn '${n.min_first_turn}' must be numeric (or blank for the package default 0.5)"
                } else if (m == 'tmm' && n.bin_size && !isNum(n.bin_size)) {
                    errs << "norm_run ${n.norm_run}: bin_size '${n.bin_size}' must be numeric (or empty for tmm_bin_size)"
                } else if (textOf(n.min_first_turn)) {
                    log.warn "min_first_turn is only used by chipseqspikeinfree; ignored"
                }
            }
        }
        def dupNames = names.findAll { s -> s && names.count(s) > 1 }.unique()
        if (dupNames) errs << "duplicate norm_run values in ${params.norm_runs}: ${dupNames.join(', ')}"
        failList("${params.norm_runs}: problems", errs)

        def runByName = nrows.collectEntries { n -> [(n.norm_run): n] }
        // IgG controls may leave norm_run empty: their bigWigs are then unscaled (multiplier 1, same units
        // as the targets). A global factor for IgG has no meaning (no enrichment to anchor on).
        def noRun = activeRows.findAll { row -> !row.norm_run && !(row.sample_id in controlIds) }.collect { row -> row.sample_id }
        if (noRun) errs << "target samples with an empty norm_run (needed unless normalisation is cpm; only IgG controls may leave it empty): ${noRun.join(', ')}"
        def unknownRun = activeRows.findAll { row -> row.norm_run && !runByName.containsKey(row.norm_run) }.collect { row -> "${row.sample_id} -> ${row.norm_run}" }
        if (unknownRun) errs << "norm_run not defined in ${params.norm_runs}: ${unknownRun.join(', ')}"
        failList("${params.samplesheet}: norm_run problems", errs)

        activeRows.each { row -> normRunOf[row.sample_id] = row.norm_run }
        nrows.each { n ->
            def members = activeRows.findAll { row -> row.norm_run == n.norm_run }.collect { row -> row.sample_id }
            def m = n.method.toLowerCase()
            if (!members) {
                log.warn "norm_run ${n.norm_run} has no active samples: dropped"
            } else {
                def bamMembers = members.findAll { s -> bamById.containsKey(s) }
                if (m == 'spikein' && bamMembers) errs << "norm_run ${n.norm_run} (spikein) contains pre-aligned BAM samples, which have no spike-in reads: ${bamMembers.join(', ')}"
                if (members.size() < 2 && m != 'cpm') errs << "norm_run ${n.norm_run} has ${members.size()} library: a normalisation run needs at least 2 libraries; use method cpm or merge runs"
                members.findAll { s -> s in controlIds }.each { s -> log.warn "control ${s} is in norm_run ${n.norm_run}; its factor has no global meaning" }
                def refCond = textOf(n.ref_condition)
                def runConds = members.collect { s -> meta[s].condition.toString() }.unique()
                if (refCond && !(refCond in runConds)) errs << "norm_run ${n.norm_run}: ref_condition '${refCond}' is not a condition of the run's samples (${runConds.join(', ')})"
                // Methods of the run: primary + norm_methods + spike-in when toggled (a 1-library run only has its cpm primary).
                def runMethods = members.size() < 2 ? [m] : ([m] + normMethods + (isTrue(params.norm_spikein) ? ['spikein'] : [])).unique()
                def unavailable = [:]   // method -> reason, for methods that are not launched
                if (('spikein' in runMethods) && m != 'spikein' && bamMembers) {
                    unavailable['spikein'] = 'pre-aligned BAM samples have no spike-in reads'
                    log.warn "norm_run ${n.norm_run}: norm_spikein is true but ${bamMembers.join(', ')} are pre-aligned BAMs (no spike-in reads): spike-in is unavailable for this run"
                }
                // ChIPseqSpikeInFree settings: the run's cells when CSIF is its primary; as a diagnostic, numeric cells else the csif_* params.
                def csif = [bin_size: n.bin_size, cutoff: n.cutoff, max_turns: n.max_turns, min_first_turn: textOf(n.min_first_turn)]
                if (m != 'chipseqspikeinfree') {
                    csif = [bin_size: (m != 'tmm' && isNum(n.bin_size)) ? textOf(n.bin_size) : textOf(params.csif_bin_size),
                            cutoff: isNum(n.cutoff) ? textOf(n.cutoff) : textOf(params.csif_cutoff),
                            max_turns: isNum(n.max_turns) ? textOf(n.max_turns) : textOf(params.csif_max_turns),
                            min_first_turn: isNum(n.min_first_turn) ? textOf(n.min_first_turn) : '']
                }
                normRuns << [norm_run: n.norm_run, method: m,
                             bin_size: (m == 'tmm' && !n.bin_size) ? params.tmm_bin_size.toString() : n.bin_size,
                             cutoff: n.cutoff, max_turns: n.max_turns, min_first_turn: textOf(n.min_first_turn), sample_ids: members,
                             methods: runMethods, unavailable: unavailable, launched: runMethods.findAll { x -> !unavailable.containsKey(x) },
                             csif: csif, ref_condition: refCond, ref_source: refCond ? 'norm_runs.csv ref_condition' : '']
            }
        }
        failList("${params.norm_runs}: problems", errs)
    }
    def normRunsById = normRuns.collectEntries { n -> [(n.norm_run): n] }

    // ---- Peak sets (csaw) -----------------------------------------------------
    def peaksets = []
    def skippedPeaksets = []
    if (peakOn('run_csaw') || peakOn('run_csaw_cutoff_scan')) {
        if (!params.peaksets) fail("peaksets is required when run_csaw or run_csaw_cutoff_scan is true (set peaksets in your config.yaml, or turn them off).")
        def prows = readTable('peaksets', params.peaksets, ['peakset', 'target_group', 'control_group', 'small_filt', 'large_filt', 'merge_bp', 'max_width', 'fdr', 'lfc'])
        if (!prows) log.warn "${params.peaksets} has no rows: no csaw peak sets will be called"
        // A peak set whose target or control group has no active sample (e.g. samples not in the samplesheet) is skipped.
        prows = prows.findAll { p ->
            def gone = ['target_group', 'control_group'].findAll { col -> p[col] && !groups.containsKey(p[col]) }.collect { col -> "${col} '${p[col]}'" }
            if (gone) {
                log.warn "${params.peaksets}: peakset ${p.peakset} skipped: ${gone.join(', ')} has no active samples in ${params.samplesheet}"
                skippedPeaksets << p.peakset
            }
            return !gone
        }
        def pnames = prows.collect { p -> p.peakset }
        def dupPeak = pnames.findAll { s -> s && pnames.count(s) > 1 }.unique()
        if (dupPeak) errs << "duplicate peakset values in ${params.peaksets}: ${dupPeak.join(', ')}"
        errs.addAll(csawDefaultsProblems())
        prows.each { p ->
            if (!p.peakset || !(p.peakset ==~ idRe)) errs << "peakset '${p.peakset}' must match ^${idRe}\$ (it is used in file names)"
            ['target_group', 'control_group'].each { col ->
                if (!groups.containsKey(p[col])) errs << "peakset ${p.peakset}: ${col} '${p[col]}' is not a group of any active (non-excluded, non-skipped) sample. Groups: ${groups.keySet().join(', ')}"
            }
            def tg = groups[p.target_group]
            def cg = groups[p.control_group]
            if (tg && cg) {
                p.target_ids = tg.findAll { s -> !(s in controlIds) }
                p.control_ids = cg.findAll { s -> s in controlIds }
                def badControls = cg.findAll { s -> !(s in controlIds) }
                if (!p.target_ids) errs << "peakset ${p.peakset}: target_group ${p.target_group} has no non-control (non-IgG) samples"
                if (!p.control_ids) errs << "peakset ${p.peakset}: control_group ${p.control_group} has no control (IgG) samples"
                if (badControls) errs << "peakset ${p.peakset}: control_group ${p.control_group} contains non-control samples: ${badControls.join(', ')}"
            }
            // Antibody of the peak set: the one antibody of its target samples (keys csaw_defaults).
            def antibodies = (p.target_ids ?: []).collect { s -> meta[s].antibody }.unique()
            if (antibodies.size() > 1) errs << "peakset ${p.peakset}: target_group ${p.target_group} has more than one antibody (${antibodies.join(', ')}); a peak set needs one antibody"
            p.antibody = antibodies.size() == 1 ? antibodies[0] : ''
            // Thresholds: non-blank row cell > csaw_defaults[antibody][col] > csaw_<col>. Each must be numeric afterwards.
            p.sources = [:]
            csawThresholdCols().each { col ->
                def v = textOf(p[col])
                def src = 'row'
                if (!v) {
                    v = textOf(csawDefaultsFor(p.antibody)[col])
                    src = 'antibody'
                }
                if (!v) {
                    v = textOf(params["csaw_${col}".toString()])
                    src = 'global'
                }
                if (!isNum(v)) errs << "peakset ${p.peakset}: ${col} '${v}' (source: ${src}) must be numeric"
                p[col] = v
                p.sources[col] = src
            }
        }
        failList("${params.peaksets}: problems", errs)
        prows.each { p ->
            log.info "csaw peak set ${p.peakset} (antibody ${p.antibody ?: 'none'}): " +
                csawThresholdCols().collect { c -> "${c}=${p[c]} (${p.sources[c]})" }.join(', ')
        }
        writePeaksetsResolved(prows)
        peaksets = prows.collect { p ->
            [peakset: p.peakset, target_group: p.target_group, control_group: p.control_group, small_filt: p.small_filt,
             large_filt: p.large_filt, merge_bp: p.merge_bp, max_width: p.max_width, fdr: p.fdr, lfc: p.lfc,
             target_ids: p.target_ids, control_ids: p.control_ids]
        }
    }

    // ---- Reference files ------------------------------------------------------
    def anyFastq = fqById.size() > 0
    def indexFiles = []
    def indexBase = null
    if (anyFastq) {
        if (!params.genome_index) {
            errs << "genome_index (bowtie2 index prefix) is required for FASTQ samples"
        } else {
            indexBase = file(hostPathNote('genome_index', params.genome_index)).name
            indexFiles = files("${hostPath(params.genome_index)}*.bt2*").findAll { f -> f.name.startsWith("${indexBase}.".toString()) }.sort { f -> f.name }
            def names = indexFiles.collect { f -> f.name }
            def missingParts = ['1', '2', '3', '4', 'rev.1', 'rev.2'].findAll { x -> !(names.contains("${indexBase}.${x}.bt2".toString()) || names.contains("${indexBase}.${x}.bt2l".toString())) }
            if (!indexFiles) errs << "genome_index ${params.genome_index}: no ${indexBase}*.bt2 / .bt2l files found"
            else if (missingParts) errs << "genome_index ${params.genome_index}: incomplete bowtie2 index, missing ${missingParts.collect { x -> indexBase + '.' + x + '.bt2(l)' }.join(', ')}"
        }
        if (!params.blacklist) errs << "blacklist (BED) is required for FASTQ samples"
    }
    ['chrom_sizes', 'blacklist', 'gtf', 'tss_bed'].each { k ->
        def required = k == 'chrom_sizes' || (k == 'tss_bed' && !params.gtf) || (k == 'blacklist' && anyFastq)
        if (!params[k]) {
            if (k == 'tss_bed' && required) errs << "tss_bed is required when gtf is not set: set gtf (preferred) or tss_bed"
            else if (required && !(k == 'blacklist' && anyFastq)) errs << "${k} is required"   // blacklist already reported above
        } else if (!file(hostPathNote(k, params[k])).exists()) {
            errs << "${k} not found: ${params[k]}"
        }
    }
    failList("reference file problems", errs)
    def refs = [index_files: indexFiles, index_base: indexBase,
                chrom_sizes: file(hostPath(params.chrom_sizes)),
                blacklist: params.blacklist ? file(hostPath(params.blacklist)) : null,
                gtf: params.gtf ? file(hostPath(params.gtf)) : null,
                tss_bed: params.tss_bed ? file(hostPath(params.tss_bed)) : null]

    // ---- Optional steps: heatmaps, tracks, annotation, intersect ---------------
    def csawNames = peakOn('run_csaw') ? peaksets.collect { p -> p.peakset } : []

    // Region classifications (run_classify): each id gives the regions <id>.retained, <id>.lost, <id>.unclassified.
    def classifications = []
    def classNames = []
    def skippedClasses = []
    def skippedIntersects = []
    if (isTrue(params.run_classify)) {
        if (!params.classifications) fail("run_classify is true but classifications (the classifications CSV path) is not set in config.yaml.")
        def crows = readTable('classifications', params.classifications, ['classification', 'regions', 'ref_group', 'treat_group', 'norm_run', 'lfc_retained', 'lfc_lost', 'fdr'])
        if (!crows) fail("${params.classifications} has no data rows.")
        // A classification whose peak set was skipped, or whose ref or treat group has no active sample, is skipped.
        crows = crows.findAll { c ->
            def why = []
            if (c.regions in skippedPeaksets) why << "regions '${c.regions}' is a skipped peak set"
            ['ref_group', 'treat_group'].each { col -> if (c[col] && !groups.containsKey(c[col])) why << "${col} '${c[col]}' has no active samples in ${params.samplesheet}" }
            if (why) {
                log.warn "${params.classifications}: classification ${c.classification} skipped: ${why.join('; ')}"
                skippedClasses.addAll(['retained', 'lost', 'unclassified'].collect { k -> "${c.classification}.${k}".toString() })
            }
            return !why
        }
        def cids = crows.collect { c -> c.classification }
        crows.each { c ->
            def lbl = "classification ${c.classification}"
            if (!(c.classification ==~ idRe)) errs << "classification '${c.classification}' must match ^${idRe}\$ (it is used in file names)"
            // Regions: a csaw peak set of this run or a BED file (classification names are not accepted here).
            def reg = classifyBed(c.regions, csawNames, [], [])
            if (!reg) errs << "${lbl}: regions '${c.regions}' is neither a csaw peak set called in this run (${regionDesc(csawNames, [], [])}) nor an existing BED file"
            // Groups: must exist and hold targets only (no IgG), one or more samples each.
            ['ref_group', 'treat_group'].each { col ->
                if (!groups.containsKey(c[col])) errs << "${lbl}: ${col} '${c[col]}' is not a group of any active (non-excluded, non-skipped) sample. Groups: ${groups.keySet().join(', ')}"
            }
            if (c.ref_group && c.ref_group == c.treat_group) errs << "${lbl}: ref_group and treat_group are both '${c.ref_group}'; they must differ"
            def refIds = groups[c.ref_group] ?: []
            def treatIds = groups[c.treat_group] ?: []
            [ref_group: refIds, treat_group: treatIds].each { col, ids ->
                if (!ids.isEmpty()) {
                    def igg = ids.findAll { s -> s in controlIds }
                    if (igg) errs << "${lbl}: ${col} '${c[col]}' contains IgG control sample(s) ${igg.join(', ')}; a classification group holds targets only"
                    if (ids.size() < 2) log.warn "${lbl}: ${col} '${c[col]}' has ${ids.size()} sample(s): no replicate dispersion: the fixed-dispersion fallback of callpeaks_csaw.R is used"
                }
            }
            // Normalisation run: the scale factors of this run are used.
            def members = (refIds + treatIds).unique()
            def normRun = textOf(c.norm_run)
            if (normRuns.isEmpty()) {
                errs << "${lbl}: needs a norm_run, but normalisation is cpm (no norm runs): the regions are scaled by the norm run's edger_norm_factor"
            } else if (!normRun) {
                def missing = members.findAll { s -> !normRunOf[s] }
                def runs = members.collect { s -> normRunOf[s] }.unique()
                if (missing) errs << "${lbl}: norm_run is blank and these samples have no norm_run: ${missing.join(', ')}"
                else if (runs.size() != 1) errs << "${lbl}: norm_run is blank and the ref and treat samples use different norm runs: " + members.collect { s -> "${s} -> ${normRunOf[s]}" }.join(', ')
                else normRun = runs[0]
            } else if (!normRunsById.containsKey(normRun)) {
                errs << "${lbl}: norm_run '${normRun}' is not defined in ${params.norm_runs}"
            } else {
                def other = members.findAll { s -> normRunOf[s] != normRun }
                if (other) errs << "${lbl}: norm_run ${normRun} is not the norm_run of these samples: " + other.collect { s -> "${s} -> ${normRunOf[s]}" }.join(', ')
            }
            // Thresholds: blank cell -> classify_* param. Numeric; fdr in (0, 1].
            def lfcRetained = textOf(c.lfc_retained) ?: textOf(params.classify_lfc_retained)
            def lfcLost = textOf(c.lfc_lost) ?: textOf(params.classify_lfc_lost)
            def fdr = textOf(c.fdr) ?: textOf(params.classify_fdr)
            [lfc_retained: lfcRetained, lfc_lost: lfcLost, fdr: fdr].each { col, v ->
                if (!isNum(v)) errs << "${lbl}: ${col} '${v}' must be numeric"
            }
            if (isNum(fdr) && !(fdr.toDouble() > 0 && fdr.toDouble() <= 1)) errs << "${lbl}: fdr '${fdr}' must be in (0, 1]"
            if (isNum(lfcRetained) && isNum(lfcLost) && lfcLost.toDouble() > lfcRetained.toDouble()) {
                log.warn "${lbl}: lfc_lost (${lfcLost}) is above lfc_retained (${lfcRetained}): check the thresholds"
            }
            classifications << [classification: c.classification, regions: reg, ref_group: c.ref_group, treat_group: c.treat_group,
                                ref_ids: refIds, treat_ids: treatIds, norm_run: normRun, lfc_retained: lfcRetained,
                                lfc_lost: lfcLost, fdr: fdr]
            classNames.addAll(['retained', 'lost', 'unclassified'].collect { k -> "${c.classification}.${k}".toString() })
        }
        def cdups = cids.findAll { s -> s && cids.count(s) > 1 }.unique()
        if (cdups) errs << "duplicate classification values in ${params.classifications}: ${cdups.join(', ')}"
        // A class region name must not be a csaw peak set name (csaw names of this run, or any peakset row).
        def peakNames = csawNames + peaksets.collect { p -> p.peakset }
        classNames.findAll { n -> n in peakNames }.each { n ->
            errs << "classification region '${n}' has the same name as a csaw peak set: rename the classification or the peak set"
        }
        failList("${params.classifications}: problems", errs)
    }

    // MA reference condition of every run with >= 2 libraries: norm_runs.csv ref_condition, else the one condition of the
    // ref_group(s) of the classifications that use the run, else the condition of the run's first library.
    normRuns.findAll { n -> n.sample_ids.size() >= 2 }.each { n ->
        if (!n.ref_condition) {
            def clsConds = classifications.findAll { c -> c.norm_run == n.norm_run }
                .collectMany { c -> c.ref_ids.collect { s -> meta[s].condition.toString() } }.unique()
            if (clsConds.size() == 1) {
                n.ref_condition = clsConds[0]
                n.ref_source = 'condition of the ref_group of the classifications using this run'
            } else {
                n.ref_condition = meta[n.sample_ids[0]].condition.toString()
                n.ref_source = clsConds ? "first library of the run: the classifications' ref_groups have several conditions (${clsConds.join(', ')})" : 'first library of the run'
            }
        }
        log.warn "norm_run ${n.norm_run}: MA reference condition = ${n.ref_condition} (${n.ref_source})"
    }
    def colorRe = '(#[0-9A-Fa-f]{6}|[A-Za-z]+)'
    def heatmaps = []
    def tracks = []
    def intersects = []

    // Intersections: records of bed_a that overlap bed_b. Parsed before the heatmaps, tracks and annotation, which may name an intersect.
    if (isTrue(params.run_intersect)) {
        if (!params.intersects) fail("run_intersect is true but intersects (the intersects CSV path) is not set in config.yaml.")
        def irows = readTable('intersects', params.intersects, ['name', 'bed_a', 'bed_b', 'label_a', 'label_b', 'min_overlap', 'reciprocal'])
        // An intersect that names a skipped peak set or classification region is skipped.
        irows = irows.findAll { r ->
            def bad = [r.bed_a, r.bed_b].findAll { b -> (b in skippedPeaksets) || (b in skippedClasses) }
            if (bad) {
                log.warn "${params.intersects}: intersect ${r.name} skipped: ${bad.join(', ')} was skipped (no active samples)"
                skippedIntersects << r.name
            }
            return !bad
        }
        def inames = irows.collect { r -> r.name }
        irows.each { r ->
            def lbl = "intersect ${r.name}"
            if (!(r.name ==~ idRe)) errs << "intersect '${r.name}' must match ^${idRe}\$ (it is used in file names)"
            // bed_a / bed_b: csaw peak sets, classification regions or BED files (not other intersect names).
            def bedA = classifyBed(r.bed_a, csawNames, classNames, [])
            if (!bedA) errs << "${lbl}: bed_a '${r.bed_a}' is not a csaw peak set called in this run, a classification region, or an existing BED file. Valid names: ${regionDesc(csawNames, classNames, [])}"
            def bedB = classifyBed(r.bed_b, csawNames, classNames, [])
            if (!bedB) errs << "${lbl}: bed_b '${r.bed_b}' is not a csaw peak set called in this run, a classification region, or an existing BED file. Valid names: ${regionDesc(csawNames, classNames, [])}"
            def minOv = r.min_overlap
            if (minOv && !(isNum(minOv) && minOv.toDouble() > 0 && minOv.toDouble() <= 1)) errs << "${lbl}: min_overlap '${minOv}' must be a number between 0 and 1 (or empty)"
            def recipRaw = r.reciprocal.toLowerCase()
            if (!(recipRaw in ['', 'true', 'false', '1', '0', 'yes', 'no'])) errs << "${lbl}: reciprocal '${r.reciprocal}' must be true or false (or empty)"
            def reciprocal = isTrue(r.reciprocal) ? 'true' : 'false'
            if (reciprocal == 'true' && !minOv) errs << "${lbl}: reciprocal needs min_overlap"
            def outputRaw = (r.output ?: '').toLowerCase()
            if (!(outputRaw in ['', 'records', 'overlap'])) errs << "${lbl}: output '${r.output}' must be records or overlap (or empty = records)"
            def outputMode = outputRaw ?: 'records'
            def labelA = r.label_a ?: bedLabel(r.bed_a, bedA)
            def labelB = r.label_b ?: bedLabel(r.bed_b, bedB)
            intersects << [name: r.name, bed_a: bedA, bed_b: bedB, label_a: labelA, label_b: labelB, min_overlap: minOv, reciprocal: reciprocal, output: outputMode]
        }
        def idups = inames.findAll { s -> s && inames.count(s) > 1 }.unique()
        if (idups) errs << "duplicate intersect values in ${params.intersects}: ${idups.join(', ')}"
        failList("${params.intersects}: problems", errs)
    }
    // Intersect names are heatmap regions, track beds and annotation extras: never the name of a peak set or classification region.
    def intersectNames = intersects.collect { r -> r.name }
    intersectNames.findAll { n -> (n in csawNames) || (n in peaksets.collect { p -> p.peakset }) || (n in classNames) }.each { n ->
        errs << "intersect '${n}' has the same name as a csaw peak set / classification region: rename it"
    }
    failList("${params.intersects}: problems", errs)

    // Heatmaps: one matrix, one plot and optional clustering per row.
    if (isTrue(params.run_heatmaps)) {
        if (!params.heatmaps) fail("run_heatmaps is true but heatmaps (the heatmaps CSV path) is not set in config.yaml.")
        // Styling params: one value for every heatmap (the deepTools engine uses them; EnrichedHeatmap does not).
        def styleErrs = []
        if (!(textOf(params.heatmap_reference_point) in ['center', 'TSS', 'TES'])) styleErrs << "heatmap_reference_point '${textOf(params.heatmap_reference_point)}' must be center, TSS or TES"
        if (!(textOf(params.heatmap_missing_as_zero).toLowerCase() in ['true', 'false', '1', '0', 'yes', 'no'])) styleErrs << "heatmap_missing_as_zero '${textOf(params.heatmap_missing_as_zero)}' must be true or false"
        ['heatmap_average_type_bins', 'heatmap_average_type_summary'].each { k ->
            if (!(textOf(params[k]) in ['mean', 'median', 'min', 'max', 'std', 'sum'])) styleErrs << "${k} '${textOf(params[k])}' must be one of mean, median, min, max, std, sum"
        }
        if (!(textOf(params.heatmap_plot_type) in ['lines', 'fill', 'se', 'std'])) styleErrs << "heatmap_plot_type '${textOf(params.heatmap_plot_type)}' must be one of lines, fill, se, std (plotHeatmap 3.5.6 --plotType)"
        def legendOk = ['best', 'upper-right', 'upper-left', 'upper-center', 'lower-left', 'lower-center', 'lower-right', 'center', 'center-left', 'center-right', 'none']
        if (!(textOf(params.heatmap_legend_location) in legendOk)) styleErrs << "heatmap_legend_location '${textOf(params.heatmap_legend_location)}' must be one of ${legendOk.join(', ')}"
        ['heatmap_height', 'heatmap_width'].each { k ->
            if (!(isNum(params[k]) && params[k].toString().toDouble() > 0)) styleErrs << "${k} '${params[k]}' must be a positive number"
        }
        ['heatmap_x_label', 'heatmap_ref_label', 'heatmap_y_label'].each { k ->
            if (hasUnsafeChar(textOf(params[k]))) styleErrs << "${k} must not contain quotes, backticks or a dollar sign"
        }
        failList("heatmap styling params", styleErrs)

        def hrows = readTable('heatmaps', params.heatmaps, ['heatmap', 'regions', 'samples', 'labels', 'window_bp', 'bin_bp', 'sort_using', 'cluster_k', 'z_min', 'z_max', 'y_max', 'colors', 'engine'])
        if (!hrows) fail("${params.heatmaps} has no data rows.")
        def skippedRegions = skippedPeaksets + skippedClasses + skippedIntersects
        hrows = hrows.collect { h -> pruneRow(params.heatmaps.toString(), params.samplesheet.toString(), 'heatmap', h, activeIds, skippedRegions,
                                              ['labels', 'colors', 'z_min', 'z_max', 'y_max'], 'regions', ['region_labels'], true) }.findAll { h -> h != null }
        def hids = hrows.collect { h -> h.heatmap }
        hrows.each { h ->
            def lbl = "heatmap ${h.heatmap}"
            if (!(h.heatmap ==~ idRe)) errs << "heatmap '${h.heatmap}' must match ^${idRe}\$ (it is used in file names)"
            // Regions: csaw peak sets, classification regions or BED files. Several share one matrix (deepTools).
            def regionNames = csvList(h.regions)
            if (!regionNames) errs << "${lbl}: regions is empty"
            def regionMaps = []
            regionNames.each { r ->
                def bm = classifyBed(r, csawNames, classNames, intersectNames)
                if (bm) regionMaps << bm
                else errs << "${lbl}: regions entry '${r}' is not a csaw peak set called in this run, a classification region, or an existing BED file. Valid names: ${regionDesc(csawNames, classNames, intersectNames)}"
            }
            def dupRegions = regionNames.findAll { s -> regionNames.count(s) > 1 }.unique()
            if (dupRegions) errs << "${lbl}: regions lists a region twice: ${dupRegions.join(', ')}"
            def stageNames = regionMaps.collect { bm -> bedStageName(bm) }
            def dupStage = stageNames.findAll { s -> stageNames.count(s) > 1 }.unique()
            if (dupStage) errs << "${lbl}: regions have the same file name (${dupStage.join(', ')}); the matrix needs a unique file name per region"
            def nR = regionNames.size()
            // region_labels: one per region; default = the region names when there are several, none for one.
            def regionLbls = csvListKeep(h.region_labels)
            if (regionLbls && regionLbls.size() != nR) errs << "${lbl}: region_labels has ${regionLbls.size()} entries but regions has ${nR}"
            regionLbls.findAll { s -> hasUnsafeChar(s) }.each { s -> errs << "${lbl}: region_labels entry '${s}' must not contain quotes, backticks or a dollar sign" }
            if (!regionLbls) regionLbls = nR > 1 ? regionNames : []
            def samples = csvList(h.samples)
            errs.addAll(sampleProblems(lbl, samples, activeIds, params.samplesheet.toString(), params.heatmaps.toString(), meta.keySet().toList()))
            def nS = samples.size()
            def labels = csvList(h.labels)
            if (!labels) labels = samples
            errs.addAll(countProblems(lbl, 'labels', labels, nS))
            def colors = csvList(h.colors)
            if (!colors) colors = samples.collect { s -> '#8c96c6' }
            errs.addAll(countProblems(lbl, 'colors', colors, nS))
            errs.addAll(colourProblems(lbl, colors, colorRe))
            def winBp = h.window_bp
            def binBp = h.bin_bp
            def winOk = (winBp ==~ /[0-9]+/) && winBp.toInteger() > 0
            def binOk = (binBp ==~ /[0-9]+/) && binBp.toInteger() > 0
            if (!winOk) errs << "${lbl}: window_bp '${winBp}' must be a positive integer"
            if (!binOk) errs << "${lbl}: bin_bp '${binBp}' must be a positive integer"
            if (winOk && binOk && winBp.toInteger() % binBp.toInteger() != 0) errs << "${lbl}: window_bp (${winBp}) must be a multiple of bin_bp (${binBp})"
            def ck = h.cluster_k
            def kParts = (ck ==~ /[0-9]+-[0-9]+/) ? ck.split('-') : []
            def clusterK = ''
            if (ck == '' || ck == '0') clusterK = ''
            else if (ck ==~ /[0-9]+/ && ck.toInteger() >= 2) clusterK = ck
            else if (kParts && kParts[0].toInteger() >= 2 && kParts[0].toInteger() < kParts[1].toInteger()) clusterK = ck
            else errs << "${lbl}: cluster_k '${ck}' must be empty, an integer >= 2, or a range like 2-6 with min < max"
            def eng = h.engine ? h.engine.toLowerCase() : 'deeptools'
            if (!(eng in ['deeptools', 'enrichedheatmap'])) errs << "${lbl}: engine '${h.engine}' must be deeptools or enrichedheatmap"
            // One region set per matrix: EnrichedHeatmap and clustering take one.
            if (nR > 1 && eng == 'enrichedheatmap') errs << "${lbl}: engine enrichedheatmap takes one region set; this row lists ${nR}"
            if (nR > 1 && clusterK) errs << "${lbl}: cluster_k clusters one region set; this row lists ${nR}"
            def zMin = csvList(h.z_min)
            def zMax = csvList(h.z_max)
            def yMax = csvList(h.y_max)
            errs.addAll(numListProblems(lbl, 'z_min', zMin, nS))
            errs.addAll(limitListProblems(lbl, 'z_max', zMax, nS, true))
            errs.addAll(limitListProblems(lbl, 'y_max', yMax, nS, true))
            // 'auto' limits are computed from the plotted matrix (HEATMAP_LIMITS): deepTools only.
            def autoLimits = (zMax + yMax).contains('auto')
            if (autoLimits && eng == 'enrichedheatmap') errs << "${lbl}: auto limits need the deepTools engine (z_max / y_max auto with engine enrichedheatmap)"
            def sortUsing = h.sort_using
            if (sortUsing && !(sortUsing in samples)) errs << "${lbl}: sort_using '${sortUsing}' is not one of its samples"
            def sortIdx = (sortUsing && (sortUsing in samples)) ? (samples.indexOf(sortUsing) + 1).toString() : ''
            def clusterSamples = csvList(h.cluster_samples)
            clusterSamples.findAll { s -> !(s in samples) }.each { s -> errs << "${lbl}: cluster_samples contains '${s}', which is not in samples" }
            if (!clusterSamples) clusterSamples = samples
            def clusterLabels = csvList(h.cluster_labels)
            def legendTitle = h.legend_title ?: ''
            def legendTitles = legendTitle ? [legendTitle] * nS : samples.collect { s -> meta[s]?.antibody }
            def antibodies = samples.collect { s -> meta[s]?.antibody }
            heatmaps << [heatmap: h.heatmap, regions: regionMaps, region_labels: regionLbls, auto: autoLimits, antibodies: antibodies,
                         samples: samples, labels: labels, colors: colors, window_bp: winBp, bin_bp: binBp, sort_using: sortUsing, sort_idx: sortIdx,
                         cluster_k: clusterK, cluster_samples: clusterSamples, cluster_labels: clusterLabels,
                         z_min: zMin, z_max: zMax, y_max: yMax, engine: eng, legend_titles: legendTitles]
        }
        def hdups = hids.findAll { s -> s && hids.count(s) > 1 }.unique()
        if (hdups) errs << "duplicate heatmap values in ${params.heatmaps}: ${hdups.join(', ')}"
        failList("${params.heatmaps}: problems", errs)
    }

    // Tracks: pyGenomeTracks figures over gene or region loci.
    if (isTrue(params.run_tracks)) {
        if (!params.tracks && !peakOn('run_peak_qc')) fail("run_tracks is true but tracks (the tracks CSV path) is not set in config.yaml (or set run_peak_qc: true for the automatic peak QC tracks only).")
        if (!params.gtf) errs << "run_tracks needs gtf (set gtf in config.yaml): gene names are looked up in it"
        // Styling params: one value for every track set. dflt = the max_value of a blank entry.
        def styleErrs = []
        if (!(textOf(params.track_number_of_bins) ==~ /[0-9]+/ && textOf(params.track_number_of_bins).toInteger() > 0)) styleErrs << "track_number_of_bins '${params.track_number_of_bins}' must be a positive integer"
        ['track_height', 'track_fontsize', 'track_bed_height', 'track_scalebar_fontsize', 'track_genes_height'].each { k ->
            if (!(isNum(params[k]) && params[k].toString().toDouble() > 0)) styleErrs << "${k} '${params[k]}' must be a positive number"
        }
        if (!(isNum(params.track_log2_min) && isNum(params.track_log2_max) && params.track_log2_min.toString().toDouble() < params.track_log2_max.toString().toDouble())) {
            styleErrs << "track_log2_min (${params.track_log2_min}) and track_log2_max (${params.track_log2_max}) must be numbers with track_log2_min < track_log2_max"
        }
        if (!(textOf(params.track_bed_color) ==~ colorRe)) styleErrs << "track_bed_color '${params.track_bed_color}' is not a #RRGGBB hex or a colour name"
        def dflt = textOf(params.tracks_max_value_default)
        if (!(dflt in ['shared', 'auto'] || isNum(dflt))) styleErrs << "tracks_max_value_default '${dflt}' must be shared, auto or a number"
        failList("track styling params", styleErrs)

        def trows = params.tracks ? readTable('tracks', params.tracks, ['track_set', 'locus_or_gene', 'flank', 'samples', 'labels', 'colors', 'max_value', 'beds', 'scalebar']) : []
        def skippedBeds = skippedPeaksets + skippedClasses + skippedIntersects
        trows = trows.collect { t -> pruneRow(params.tracks.toString(), params.samplesheet.toString(), 'track_set', t, activeIds, skippedBeds,
                                              ['labels', 'colors', 'max_value'], 'beds', ['bed_labels', 'bed_colors'], false) }.findAll { t -> t != null }
        def tids = trows.collect { t -> t.track_set }
        def geneRe = '^[A-Za-z0-9][A-Za-z0-9_.-]*$'
        def regionRe = '^[^:\\s,]+:[0-9,]+-[0-9,]+$'
        trows.each { t ->
            def lbl = "track_set ${t.track_set}"
            if (!(t.track_set ==~ idRe)) errs << "track_set '${t.track_set}' must match ^${idRe}\$ (it is used in file names)"
            def loci = csvList(t.locus_or_gene)
            if (!loci) errs << "${lbl}: locus_or_gene is empty"
            loci.findAll { x -> !(x ==~ geneRe) && !(x ==~ regionRe) }.each { x -> errs << "${lbl}: locus '${x}' is neither a gene name (letters, digits, _ . -) nor a region chr:start-end" }
            def geneLoci = loci.findAll { x -> x ==~ geneRe }
            def flanks = csvList(t.flank)
            def flankOut = []
            if (flanks && !(flanks.size() in [1, loci.size()])) errs << "${lbl}: flank has ${flanks.size()} entries; give one value or one per locus (${loci.size()})"
            else if (!flanks && geneLoci) errs << "${lbl}: flank is empty but locus_or_gene has a gene name"
            else {
                flanks.findAll { f -> !(f ==~ /[0-9]+/) }.each { f -> errs << "${lbl}: flank '${f}' must be a non-negative integer" }
                flankOut = flanks ? (flanks.size() == 1 ? loci.collect { x -> flanks[0] } : flanks) : loci.collect { x -> '0' }
            }
            def samples = csvList(t.samples)
            errs.addAll(sampleProblems(lbl, samples, activeIds, params.samplesheet.toString(), params.tracks.toString(), meta.keySet().toList()))
            def nS = samples.size()
            def labels = csvList(t.labels)
            if (!labels) labels = samples
            errs.addAll(countProblems(lbl, 'labels', labels, nS))
            def colors = csvList(t.colors)
            if (!colors) colors = samples.collect { s -> '#8c96c6' }
            errs.addAll(countProblems(lbl, 'colors', colors, nS))
            errs.addAll(colourProblems(lbl, colors, colorRe))
            // max_value: none, one, or one per sample; each blank, auto, shared, number or group:NAME.
            def maxEntries = csvListKeep(t.max_value)
            if (!(maxEntries.size() in [0, 1, nS])) errs << "${lbl}: max_value has ${maxEntries.size()} entries; give none, one, or one per sample (${nS})"
            maxEntries.findAll { e -> !(e == '' || e == 'auto' || e == 'shared' || isNum(e) || (e.startsWith('group:') && e.substring(6) ==~ /[A-Za-z0-9_.-]+/)) }.each { e ->
                errs << "${lbl}: max_value entry '${e}' must be a number, auto, shared, group:NAME or blank"
            }
            // Beds: csaw peak sets, classification regions or BED files; optional colours and labels, one per bed.
            def bedNames = csvList(t.beds)
            def bedMaps = []
            def bedLbls = []
            def bedLabelsGiven = csvListKeep(t.bed_labels)
            if (bedLabelsGiven && bedLabelsGiven.size() != bedNames.size()) errs << "${lbl}: bed_labels has ${bedLabelsGiven.size()} entries but beds has ${bedNames.size()}"
            bedNames.eachWithIndex { b, i ->
                def bm = classifyBed(b, csawNames, classNames, intersectNames)
                if (!bm) {
                    errs << "${lbl}: beds entry '${b}' is not a csaw peak set called in this run, a classification region, or an existing BED file. Valid names: ${regionDesc(csawNames, classNames, intersectNames)}"
                    return
                }
                bedMaps << bm
                def given = bedLabelsGiven.size() == bedNames.size() ? bedLabelsGiven[i] : ''
                bedLbls << (given ?: bedLabel(b, bm))
            }
            bedLbls.findAll { s -> hasUnsafeChar(s) }.each { s -> errs << "${lbl}: bed_labels entry '${s}' must not contain quotes, backticks or a dollar sign" }
            def bedColours = csvListKeep(t.bed_colors)
            if (bedColours && bedColours.size() != bedNames.size()) errs << "${lbl}: bed_colors has ${bedColours.size()} entries but beds has ${bedNames.size()}"
            errs.addAll(colourProblems(lbl, bedColours.findAll { c -> c }, colorRe))
            def bedColourList = bedColours.size() == bedNames.size() ? bedColours.collect { c -> c ?: textOf(params.track_bed_color) } : bedNames.collect { c -> textOf(params.track_bed_color) }
            def scalebar = t.scalebar
            if (scalebar && scalebar != 'auto' && !(scalebar ==~ /[0-9]+/ && scalebar.toInteger() > 0)) errs << "${lbl}: scalebar '${scalebar}' must be empty, auto or a positive integer (bp)"
            def iniArgs = (t.ini_args ?: '').toString().trim()
            if (iniArgs && !(iniArgs.startsWith('--') && iniArgs ==~ /[A-Za-z0-9_.:=, -]+/)) errs << "${lbl}: ini_args '${iniArgs}' must start with -- and use only letters, digits and _ . : = , - (make_tracks_ini.py options, e.g. --gene-rows 2 --gap-after ID=0.5)"
            def log2ratio = isTrue(t.log2ratio)
            if (log2ratio && !isTrue(params.run_log2ratio)) errs << "${lbl}: log2ratio is true but run_log2ratio is false in config.yaml"
            def antibodies = samples.collect { s -> meta[s]?.antibody }
            def maxTokens = trackMaxTokens(maxEntries, nS, dflt, antibodies)
            tracks << [track_set: t.track_set, loci: loci, flanks: flankOut, samples: samples, labels: labels, colors: colors,
                       max_tokens: maxTokens, beds: bedMaps, bed_colors: bedColourList, bed_labels: bedLbls, scalebar: scalebar, log2ratio: log2ratio, ini_args: iniArgs]
        }
        def tdups = tids.findAll { s -> s && tids.count(s) > 1 }.unique()
        if (tdups) errs << "duplicate track_set values in ${params.tracks}: ${tdups.join(', ')}"
        if (peakOn('run_peak_qc') && ('peak_qc' in tids)) errs << "track_set 'peak_qc' is reserved for the automatic peak QC tracks (run_peak_qc): rename it in ${params.tracks}"
        failList("${params.tracks}: problems", errs)
    }

    // Peak QC (run_peak_qc): optional parameters.
    if (peakOn('run_peak_qc')) {
        def pqErrs = []
        def pqMax = textOf(params.peak_qc_max_regions)
        if (!(pqMax ==~ /[0-9]+/ && pqMax.toInteger() >= 1)) pqErrs << "peak_qc_max_regions '${pqMax}' must be a positive integer"
        ['peak_qc_loci', 'peak_qc_negative', 'peak_qc_shape_marks'].each { k ->
            if (hasUnsafeChar(textOf(params[k]))) pqErrs << "${k} must not contain quotes, backticks or a dollar sign"
        }
        if (textOf(params.peak_qc_nopeak_igg_max_fold) && !(isNum(params.peak_qc_nopeak_igg_max_fold) && params.peak_qc_nopeak_igg_max_fold.toString().toDouble() > 0)) {
            pqErrs << "peak_qc_nopeak_igg_max_fold '${params.peak_qc_nopeak_igg_max_fold}' must be a positive number (or null / '' for off)"
        }
        csvList(params.peak_qc_exclude_bed).findAll { b -> !file(hostPathNote('peak_qc_exclude_bed', b)).exists() }.each { b -> pqErrs << "peak_qc_exclude_bed: ${b} not found" }
        failList("peak QC parameters: problems", pqErrs)
    }

    // Peak annotation: needs the GTF and at least one peak source (csaw, MACS3, SEACR or annotation_extra_beds).
    def annotationExtra = []   // [label:, bm:] of annotation_extra_beds (run_annotation)
    if (peakOn('run_annotation')) {
        if (!params.gtf) errs << "run_annotation needs gtf (set gtf in config.yaml)"
        // Extra regions annotated next to the peak sets; a csaw peak set that run_csaw already annotates is skipped.
        csvList(params.annotation_extra_beds).each { ref ->
            def bm = classifyBed(ref, csawNames, classNames, intersectNames)
            if (!bm) {
                errs << "annotation_extra_beds entry '${ref}' is not a csaw peak set called in this run, a classification region, an intersect name or an existing BED file. Valid names: ${regionDesc(csawNames, classNames, intersectNames)}"
                return
            }
            if (bm.kind == 'peakset' && peakOn('run_csaw')) {
                log.warn "annotation_extra_beds: '${ref}' is a csaw peak set that is already annotated (skipped)"
                return
            }
            annotationExtra << [label: bedLabel(ref, bm), bm: bm]
        }
        // Labels name the annotation outputs: unique, also after replacing [^A-Za-z0-9_] by _ (peak_annotation.R slugs them).
        def annotLabels = (peakOn('run_csaw') ? csawNames : []) + annotationExtra.collect { x -> x.label }
        def labelDups = annotLabels.findAll { s -> annotLabels.count(s) > 1 }.unique()
        if (labelDups) errs << "annotation_extra_beds: duplicate label(s) ${labelDups.join(', ')}"
        def slugDups = annotLabels.findAll { s -> annotLabels.findAll { t -> t.replaceAll('[^A-Za-z0-9_]', '_') == s.replaceAll('[^A-Za-z0-9_]', '_') }.size() > 1 }.unique()
        if (slugDups && !labelDups) errs << "annotation labels collide after replacing characters other than A-Za-z0-9_ by _: ${slugDups.join(', ')}"
        if (!((peakOn('run_csaw') && peaksets) || runMacs3 || runSeacr || annotationExtra)) errs << "run_annotation needs at least one peak source (run_csaw with peaksets, run_macs3 or run_seacr)"
        csvList(params.annotation_colors).findAll { c -> !(c ==~ '^[^=]+=#[0-9A-Fa-f]{6}$') }.each { c -> errs << "annotation_colors: '${c}' is not label=#RRGGBB" }
        failList("annotation: problems", errs)
    }

    // ---- Summary --------------------------------------------------------------
    def on = []
    ['trim', 'remove_duplicates', 'spikein_dedup', 'run_csaw', 'run_csaw_cutoff_scan', 'run_classify', 'run_macs3', 'run_seacr', 'run_log2ratio', 'run_report',
     'run_heatmaps', 'run_tracks', 'run_peak_qc', 'run_annotation', 'run_intersect'].each { k ->
        if (k == 'run_macs3' ? runMacs3 : k == 'run_seacr' ? runSeacr : (k in ['run_csaw', 'run_csaw_cutoff_scan', 'run_peak_qc', 'run_annotation'] ? peakOn(k) : isTrue(params[k]))) on << k
    }
    def optLines = []
    if (isTrue(params.run_classify)) optLines << "classifications : ${classifications.size()} (${classifications.collect { c -> c.classification }.join(', ')})"
    if (isTrue(params.run_heatmaps)) optLines << "heatmaps       : ${heatmaps.size()} (${heatmaps.collect { h -> h.heatmap }.join(', ')})"
    if (isTrue(params.run_tracks)) optLines << "track sets     : ${tracks.size()}"
    if (isTrue(params.run_intersect)) optLines << "intersects     : ${intersects.size()}"
    def optText = optLines.collect { l -> "\n      ${l}" }.join('')
    log.info """
    CUTANDRUN  ~  ${params.project_name}
    -----------------------------------------
      libraries      : ${fqById.size()} FASTQ, ${bamById.size()} pre-aligned BAM, ${skippedIds.size()} skipped (skip_groups), ${excludedIds.size()} excluded
      groups         : ${groups.collect { g, ids -> "${g} (${ids.size()})" }.join(', ')}
      controls       : ${controlIds.size()} IgG
      assay          : ${params.experiment} (${adapterSet}), bowtie2 ${params.bowtie2_mode}
      normalisation  : ${params.normalisation}${normRuns ? ' - runs: ' + normRuns.collect { n -> "${n.norm_run} (primary ${n.method}, ${n.sample_ids.size()} samples; methods ${n.launched.join('+')})" }.join(', ') : ''}
      peak sets      : ${callPeaks() ? (peaksets ? peaksets.collect { p -> p.peakset }.join(', ') : '(none)') : '(peak calling off: call_peaks is false)'}
      steps on       : ${on.join(', ')}${optText}
      outdir         : ${params.outdir}
    """

    return [scripts: stepScripts(), images: processContainers(), adapterSet: adapterSet,
            rows: rows, activeRows: activeRows, activeIds: activeIds, skipped: skippedIds, skippedIds: skippedIds,
            excludedIds: excludedIds, meta: meta, fastqRows: fastqRows, bamRows: bamRows, fqById: fqById, bamById: bamById,
            isControl: controlIds, controlOf: controlOf, runMacs3: runMacs3, runSeacr: runSeacr, callModes: callModes.rows,
            macsBroad: callModes.macsBroad, seacrMode: callModes.seacrMode, groups: groups, groupsOf: groupsOf,
            normRuns: normRuns, normRunsById: normRunsById, normRunOf: normRunOf, peaksets: peaksets, refs: refs,
            classifications: classifications, classNames: classNames,
            heatmaps: heatmaps, tracks: tracks, intersects: intersects, annotation: peakOn('run_annotation'), annotationExtra: annotationExtra]
}

// -----------------------------------------------------------------------
// DAG helpers. A "lookup" is a value channel holding ONE map
// sample_id -> [bam, bai] (see bamLookup); per-group / per-run inputs are
// built from it with .combine(lookup).map { ... } so that the id lists, the
// staged file names and the BAM lists of one task always come from one map.
// -----------------------------------------------------------------------
def bamLookup(ch) {
    return ch.toList().map { rows -> rows.collectEntries { row -> [(row[0]): [row[1], row[2]]] } }
}

// Value channel of one map sample_id -> file, from a channel of tuple(sample_id, file).
def fileLookup(ch) {
    return ch.toList().map { rows -> rows.collectEntries { row -> [(row[0]): row[1]] } }
}

def pick(Map lk, String id, String step) {
    if (!lk.containsKey(id)) {
        throw new Exception("${step}: no file for sample '${id}' (it was not processed: check skip_groups, exclude and the samplesheet)")
    }
    return lk[id]
}

def bamsOf(Map lk, List ids, String step) { return ids.collect { id -> pick(lk, id, step)[0] } }
def baisOf(Map lk, List ids, String step) { return ids.collect { id -> pick(lk, id, step)[1] } }

def namesOf(List ids, String suffix) { return ids.collect { id -> "${id}${suffix}".toString() } }

// File of a classifyBed() map: an existing BED, or a csaw peak set / classification region (lk = one merged lookup of both).
def bedLookupInput(Map ref, Map lk, String step) {
    return ref.kind == 'path' ? file(hostPath(ref.value)) : pick(lk, ref.value, step)
}

// File name a classifyBed() map has when staged: <peakset>.peaks.bed, <classification region>.bed, or the BED's own name.
def bedStageName(Map bm) {
    if (bm.kind == 'peakset') return "${bm.value}.peaks.bed".toString()
    if (bm.kind == 'class') return "${bm.value}.bed".toString()
    if (bm.kind == 'intersect') return "${bm.value}.intersect.bed".toString()
    return file(hostPath(bm.value)).name
}

// One z_max / y_max list expanded to one value per sample; each 'auto' takes the HEATMAP_LIMITS limit of its sample
// (lims: the limits in sample order). Lists without 'auto' are returned as they are (one entry expanded to n).
def resolveAuto(List vals, int n, List lims, String what) {
    if (!vals) return []
    def full = vals.size() == 1 ? (0..<n).collect { vals[0] } : vals
    if (!full.contains('auto')) return full
    if (lims.size() < n) throw new Exception("${what}: auto limits missing: HEATMAP_LIMITS gave ${lims.size()} value(s) for ${n} sample(s)")
    return (0..<n).collect { i -> full[i] == 'auto' ? lims[i] : full[i] }
}

// Labels of the heatmap's k cluster regions: the given ones when there are k, else "cluster 1".."cluster k".
def regionLabels(List given, int k) {
    return given.size() == k ? given : (1..k).collect { i -> "cluster ${i}".toString() }
}

// Staged path of track BED i (0-based) inside TRACKS_INI's directory: tracks_bed_<i+1>/<file name>.
def stageBedName(int i, f) {
    return "tracks_bed_${i + 1}/${f.name}".toString()
}

// Cluster BEDs of one heatmap in cluster-number order (cluster1.bed, cluster2.bed, ...).
def kSorted(List beds) {
    return beds.sort(false) { b ->   // sort(false): a copy; in-place sort of a task output list races with publishing (ConcurrentModificationException)
         (b.name =~ 'cluster([0-9]+)\\.bed$')[0][1].toInteger() }
}

// Bin sizes of the BIN_COUNTS tasks (one task each, as text): tmm_bin_size when there are norm_runs (TMM diagnostics, NORM_QC),
// the bin_size of every TMM-primary run, and the replicate-concordance size (repcor_bin_size, else tmm_bin_size).
def binCountSizes(Map cfg) {
    def sizes = []
    if (cfg.normRuns) {
        sizes << textOf(params.tmm_bin_size)
        cfg.normRuns.findAll { n -> n.method == 'tmm' }.each { n -> sizes << textOf(n.bin_size) }
    }
    if (isTrue(params.run_replicate_correlation)) sizes << (textOf(params.repcor_bin_size) ?: textOf(params.tmm_bin_size))
    return sizes.unique()
}

// Input tuple of CSAW and CSAW_CUTOFF_SCAN for one peak set.
def csawInput(Map ps, Map lk, String step) {
    return [ps.peakset,
            ps.target_ids.join(','), namesOf(ps.target_ids, '.final.bam').join(','), bamsOf(lk, ps.target_ids, step), baisOf(lk, ps.target_ids, step),
            ps.control_ids.join(','), namesOf(ps.control_ids, '.final.bam').join(','), bamsOf(lk, ps.control_ids, step), baisOf(lk, ps.control_ids, step),
            ps.small_filt, ps.large_filt, ps.merge_bp, ps.max_width, ps.fdr, ps.lfc]
}

workflow {

    def cfg = prepareRun()
    reportHostPaths()

    // Each step's script: bin/<name>, or this project's replacement (see stepScripts()).
    def sc = cfg.scripts
    def theme_r = Channel.value(sc.my_theme)
    def theme_py = Channel.value(sc.my_theme_py)
    def samplesheet = Channel.value(file(hostPath(params.samplesheet)))

    // Channels that optional blocks add to (all stay empty when their block is off).
    def final_bams = Channel.empty()        // tuple(id, bam, bai), ACTIVE samples only
    def final_flagstat = Channel.empty()    // S.final.flagstat.txt
    def final_idxstats = Channel.empty()
    def spike_bams = Channel.empty()        // tuple(id, bam, bai), FASTQ samples only
    def spike_idxstats = Channel.empty()    // S.spikein.idxstats.txt, FASTQ samples only
    def filtered_bams = Channel.empty()     // tuple(id, bam, bai), pre-dedup, FASTQ samples only
    def raw_stats = Channel.empty()         // raw flagstat + idxstats
    def markdup_metrics = Channel.empty()
    def fastqc_zips = Channel.empty()
    def trim_logs = Channel.empty()
    def align_logs = Channel.empty()

    // ---- A. Preprocess ---------------------------------------------------------
    if (cfg.fastqRows) {
        def fastqc_in = Channel.fromList(cfg.fastqRows.collectMany { row -> [tuple(row[0], row[1]), tuple(row[0], row[2])] })
        FASTQC(fastqc_in, sc.fastqc)
        fastqc_zips = FASTQC.out.zip.map { sid, zip -> zip }

        def reads = Channel.fromList(cfg.fastqRows.collect { row -> tuple(row[0], row[1], row[2]) })
        if (isTrue(params.trim)) {
            TRIM(reads, cfg.adapterSet, params.trim_minlen, sc.trim)
            reads = TRIM.out.reads
            trim_logs = TRIM.out.log
        }
        ALIGN(reads, Channel.value(cfg.refs.index_files), cfg.refs.index_base, params.bowtie2_mode, (params.bowtie2_minins ?: ''), (params.bowtie2_maxins ?: ''), sc.align)
        align_logs = ALIGN.out.log
        BAM_FILTER(ALIGN.out.bam, params.mapq, sc.bam_filter)
        filtered_bams = BAM_FILTER.out.bam
        raw_stats = BAM_FILTER.out.raw_flagstat.mix(BAM_FILTER.out.raw_idxstats)
        MARKDUP(BAM_FILTER.out.bam, sc.markdup)
        markdup_metrics = MARKDUP.out.metrics

        def dedup_in = MARKDUP.out.bam.map { sid, bam, bai -> tuple(sid, bam, bai, cfg.isControl.contains(sid).toString()) }
        DEDUP_SPLIT(dedup_in, isTrue(params.remove_duplicates).toString(), isTrue(params.spikein_dedup).toString(), params.spike_prefix, params.mito_chrom, sc.dedup_split)
        spike_bams = DEDUP_SPLIT.out.spikein
        spike_idxstats = DEDUP_SPLIT.out.spike_idxstats.map { f -> tuple(f.name.replace('.spikein.idxstats.txt', ''), f) }
        BLACKLIST_FILTER(DEDUP_SPLIT.out.target, Channel.value(cfg.refs.blacklist), sc.blacklist_filter)
        BAM_FINALIZE(BLACKLIST_FILTER.out.bam, sc.bam_finalize)
        final_bams = final_bams.mix(BAM_FINALIZE.out.bam)
        final_flagstat = final_flagstat.mix(BAM_FINALIZE.out.flagstat)
        final_idxstats = final_idxstats.mix(BAM_FINALIZE.out.idxstats)
    }
    if (cfg.bamRows) {
        // one BAM_PREALIGNED call for the pre-aligned rows
        BAM_PREALIGNED(Channel.fromList(cfg.bamRows.collect { row -> tuple(row[0], row[1]) }), sc.bam_prealigned)
        final_bams = final_bams.mix(BAM_PREALIGNED.out.bam.filter { sid, bam, bai -> sid in cfg.activeIds })
        final_flagstat = final_flagstat.mix(BAM_PREALIGNED.out.flagstat.filter { f -> f.name.replace('.final.flagstat.txt', '') in cfg.activeIds })
        final_idxstats = final_idxstats.mix(BAM_PREALIGNED.out.idxstats.filter { f -> f.name.replace('.final.idxstats.txt', '') in cfg.activeIds })
    }

    // ---- B. Lookups ------------------------------------------------------------
    def bam_lookup = bamLookup(final_bams)       // id -> [bam, bai]
    def spike_idx_lookup = fileLookup(spike_idxstats)   // id -> spike-in idxstats
    def final_idx_lookup = fileLookup(final_idxstats.map { f -> tuple(f.name.replace('.final.idxstats.txt', ''), f) })   // id -> target idxstats
    def flagstat_lookup = fileLookup(final_flagstat.map { f -> tuple(f.name.replace('.final.flagstat.txt', ''), f) })

    // Tables and plots the report collects (optional blocks add to them).
    def report_tables = Channel.empty()
    def report_plots = Channel.empty()

    // ---- C. QC -----------------------------------------------------------------
    if (cfg.activeIds) {
        def align_stats_in = raw_stats.mix(final_flagstat, final_idxstats, markdup_metrics).collect()
        QC_ALIGNMENT_TABLE(align_stats_in, samplesheet, params.spike_prefix, params.mito_chrom, sc.qc_alignment_table)
        QC_ALIGNMENT_PLOT(QC_ALIGNMENT_TABLE.out.table, theme_r, sc.qc_alignment_plot)
        report_tables = report_tables.mix(QC_ALIGNMENT_TABLE.out.table)
        report_plots = report_plots.mix(QC_ALIGNMENT_PLOT.out.svg)

        def frag_in = final_bams.flatMap { row -> [row[1], row[2]] }.collect()
        QC_FRAGMENT_SIZES_TABLE(frag_in, samplesheet, sc.qc_fragment_sizes_table)
        QC_FRAGMENT_SIZES_PLOT(QC_FRAGMENT_SIZES_TABLE.out.summary, QC_FRAGMENT_SIZES_TABLE.out.distribution, theme_r, sc.qc_fragment_sizes_plot)
        report_tables = report_tables.mix(QC_FRAGMENT_SIZES_TABLE.out.summary, QC_FRAGMENT_SIZES_TABLE.out.distribution)
        report_plots = report_plots.mix(QC_FRAGMENT_SIZES_PLOT.out.svg)

        def preseq_curves = Channel.empty()
        if (cfg.fastqRows) {
            // Library complexity and preseq need the pre-dedup BAMs, so FASTQ samples only.
            LIBCOMPLEXITY_FRAGS(filtered_bams, sc.library_complexity_frags)
            PICARD_ELC(filtered_bams, sc.picard_elc)
            def libcomplexity_in = LIBCOMPLEXITY_FRAGS.out.tsv.mix(PICARD_ELC.out.metrics, markdup_metrics).collect()
            QC_LIBCOMPLEXITY_TABLE(libcomplexity_in, samplesheet, sc.qc_library_complexity_table)
            QC_LIBCOMPLEXITY_PLOT(QC_LIBCOMPLEXITY_TABLE.out.table, theme_py, sc.qc_library_complexity_plot)
            report_tables = report_tables.mix(QC_LIBCOMPLEXITY_TABLE.out.table)
            report_plots = report_plots.mix(QC_LIBCOMPLEXITY_PLOT.out.svg)

            PRESEQ(filtered_bams, sc.preseq_ccurve)
            QC_PRESEQ_PLOT(PRESEQ.out.curve.collect(), theme_r, sc.qc_preseq_plot)
            preseq_curves = PRESEQ.out.curve
            report_plots = report_plots.mix(QC_PRESEQ_PLOT.out.svg)
        }

        // One fingerprint per group (every group with at least one sample).
        def fingerprint_in = Channel.fromList(cfg.groups.collect { g, ids -> [group: g, ids: ids] })
            .combine(bam_lookup)
            .map { gr, lk -> tuple(gr.group, gr.ids.join(','), bamsOf(lk, gr.ids, 'FINGERPRINT'), baisOf(lk, gr.ids, 'FINGERPRINT')) }
        FINGERPRINT(fingerprint_in, sc.qc_fingerprint)
        report_plots = report_plots.mix(FINGERPRINT.out.svg)

        // Correlation: bams, bais and labels from ONE map, in activeIds order.
        if (isTrue(params.run_deeptools_correlation)) {
            def corr_in = bam_lookup.map { lk -> [bamsOf(lk, cfg.activeIds, 'CORRELATION'), baisOf(lk, cfg.activeIds, 'CORRELATION'), cfg.activeIds.join(',')] }
            CORRELATION(corr_in.map { c -> c[0] }, corr_in.map { c -> c[1] }, corr_in.map { c -> c[2] }, sc.qc_correlation)
            report_plots = report_plots.mix(CORRELATION.out.svg)
        }

        def multiqc_in = fastqc_zips.mix(trim_logs, align_logs, markdup_metrics, raw_stats, final_flagstat, final_idxstats, preseq_curves).collect()
        MULTIQC(multiqc_in, sc.multiqc)
    }

    // ---- C2. Bin counts and replicate concordance --------------------------------
    // One csaw pass over every active library (IgG included) per distinct bin size; TMM, NORM_QC and the replicate
    // concordance each take the columns they need. Runs also with normalisation: cpm (the replicate step needs it).
    def bc_counts = Channel.empty()   // tuple(bin_size, counts.tsv.gz, libsizes.tsv)
    def bin_sizes = binCountSizes(cfg)
    if (bin_sizes) {
        def bc_in = Channel.fromList(bin_sizes).combine(bam_lookup).map { bs, lk ->
            tuple(bs, cfg.activeIds.join(','), namesOf(cfg.activeIds, '.final.bam').join(','),
                  bamsOf(lk, cfg.activeIds, 'BIN_COUNTS'), baisOf(lk, cfg.activeIds, 'BIN_COUNTS'))
        }
        def bc_idx = final_idx_lookup.map { fk -> cfg.activeIds.collect { id -> pick(fk, id, 'BIN_COUNTS (idxstats)') } }
        BIN_COUNTS(bc_in, bc_idx, Channel.value(cfg.refs.blacklist ?: []), sc.bin_counts)
        bc_counts = BIN_COUNTS.out.counts
    }

    if (isTrue(params.run_replicate_correlation)) {
        def rc_size = textOf(params.repcor_bin_size) ?: textOf(params.tmm_bin_size)
        def rc_counts = bc_counts.filter { bs, counts, libsizes -> bs == rc_size }
        REPLICATE_CORRELATION(rc_counts, samplesheet, cfg.activeIds.join(','), sc.qc_replicate_correlation)
        REPLICATE_CORRELATION_PLOT(REPLICATE_CORRELATION.out.table, REPLICATE_CORRELATION.out.pairs_all, rc_counts, theme_r, sc.qc_replicate_correlation_plot)
        report_tables = report_tables.mix(REPLICATE_CORRELATION.out.table)
        report_plots = report_plots.mix(REPLICATE_CORRELATION_PLOT.out.svg)
    }

    // ---- D. Normalisation (skipped for normalisation: cpm) ---------------------
    // Every run computes its primary method (norm_runs.csv `method`) plus norm_methods (and spike-in with norm_spikein).
    // Only the primary feeds the bigWigs and the classification (sf_all); the others are diagnostics for NORM_QC.
    def sf_methods = Channel.empty()   // tuple(norm_run, method, <norm_run>.<method>.scale_factors.tsv): every method of every run
    def sf_all = Channel.empty()       // tuple(norm_run, scale_factors.tsv): the PRIMARY method's table only
    if (cfg.normRuns) {
        def csf_runs = cfg.normRuns.findAll { n -> 'chipseqspikeinfree' in n.launched }
        def spike_runs = cfg.normRuns.findAll { n -> 'spikein' in n.launched }
        def tmm_runs = cfg.normRuns.findAll { n -> 'tmm' in n.launched }
        def cpm_runs = cfg.normRuns.findAll { n -> 'cpm' in n.launched }

        if (csf_runs) {
            def csf_in = Channel.fromList(csf_runs).combine(bam_lookup).map { n, lk ->
                tuple(n.norm_run, n.sample_ids.join(','), namesOf(n.sample_ids, '.final.bam').join(','),
                      bamsOf(lk, n.sample_ids, 'CHIPSEQSPIKEINFREE'), baisOf(lk, n.sample_ids, 'CHIPSEQSPIKEINFREE'),
                      n.method == 'chipseqspikeinfree' ? 'primary' : 'diagnostic',
                      n.csif.bin_size, n.csif.cutoff, n.csif.max_turns, n.csif.min_first_turn)
            }
            CHIPSEQSPIKEINFREE(csf_in, Channel.value(cfg.refs.chrom_sizes), sc.norm_chipseqspikeinfree)
            sf_methods = sf_methods.mix(CHIPSEQSPIKEINFREE.out.sf)
            report_plots = report_plots.mix(CHIPSEQSPIKEINFREE.out.plot)
        }
        if (spike_runs) {
            def spike_in = Channel.fromList(spike_runs).combine(final_idx_lookup).combine(spike_idx_lookup).map { n, fk, sk ->
                tuple(n.norm_run, n.sample_ids.join(','), n.method == 'spikein' ? 'primary' : 'diagnostic',
                      n.sample_ids.collect { id -> pick(fk, id, 'SPIKEIN (target idxstats)') },
                      n.sample_ids.collect { id -> pick(sk, id, 'SPIKEIN (spike-in idxstats)') })
            }
            SPIKEIN(spike_in, params.spike_prefix, sc.norm_spikein)
            sf_methods = sf_methods.mix(SPIKEIN.out.sf)
        }
        if (tmm_runs) {
            // A TMM primary run counts at its own bin_size, a TMM diagnostic at tmm_bin_size; the counts are joined by bin size.
            def tmm_in = Channel.fromList(tmm_runs.collect { n ->
                    tuple(n.method == 'tmm' ? textOf(n.bin_size) : textOf(params.tmm_bin_size), n.norm_run, n.sample_ids.join(','), n.method == 'tmm' ? 'primary' : 'diagnostic')
                })
                .combine(bc_counts, by: 0)   // tuple(bin_size, norm_run, ids, role, counts, libsizes)
            TMM(tmm_in, sc.norm_tmm)
            sf_methods = sf_methods.mix(TMM.out.sf)
        }
        if (cpm_runs) {
            NORM_CPM(Channel.fromList(cpm_runs.collect { n -> tuple(n.norm_run, n.sample_ids.join(',')) }), sc.norm_cpm)
            sf_methods = sf_methods.mix(NORM_CPM.out.sf)
        }

        // Only the primary's table goes on to the bigWigs and the classification, so a diagnostic never delays them.
        def primary_of = cfg.normRuns.collectEntries { n -> [(n.norm_run): n.method] }
        sf_all = sf_methods.filter { nr, m, tsv -> m == primary_of[nr] }.map { nr, m, tsv -> tuple(nr, tsv) }

        // One table per run: <run>.scale_factors.tsv (byte copy of the primary's) and <run>.scale_factors_all.tsv (all methods).
        // groupKey: a run is complete when every launched method has written its table; methods that were not launched
        // (spike-in with pre-aligned BAMs) are passed as 'method=reason' and written as unavailable.
        def n_tables = cfg.normRuns.collectEntries { n -> [(n.norm_run): n.launched.size()] }
        def run_by_name = cfg.normRunsById
        def combine_in = sf_methods
            .map { nr, m, tsv -> tuple(groupKey(nr, n_tables[nr]), m, tsv) }
            .groupTuple()
            .map { k, ms, tsvs ->
                def nr = k.toString()
                def n = run_by_name[nr]
                def tables = [ms, tsvs].transpose().sort(false) { p -> p[0] }
                tuple(nr, n.method, n.methods.join(','), n.unavailable.collect { m, why -> "${m}=${why}".toString() }, tables.collect { p -> p[1] })
            }
        NORM_COMBINE(combine_in, sc.norm_combine)

        // Normalisation QC (runs with at least 2 libraries): MA pairs and reference lines from the wide table and the shared bin counts.
        def qc_runs = cfg.normRuns.findAll { n -> n.sample_ids.size() >= 2 }
        def qc_names = qc_runs.collect { n -> n.norm_run }
        if (qc_runs) {
            def qc_info = Channel.fromList(qc_runs.collect { n -> tuple(n.norm_run, n) })
            def wide_qc = NORM_COMBINE.out.wide.filter { nr, tsv -> nr in qc_names }
            def tmm_counts = bc_counts.filter { bs, counts, libsizes -> bs == textOf(params.tmm_bin_size) }
            def norm_qc_in = wide_qc
                .join(qc_info, failOnMismatch: true, failOnDuplicate: true)
                .combine(tmm_counts)
                .map { nr, wide, n, bs, counts, libsizes -> tuple(nr, n.sample_ids.join(','), n.ref_condition ?: '', n.method, wide, counts, libsizes) }
            NORM_QC(norm_qc_in, samplesheet, sc.norm_qc)
            report_tables = report_tables.mix(NORM_QC.out.tsv)

            // Normalisation MA plots: pairs table and wide table of one norm_run, joined by norm_run, with the counts.
            def norm_plot_in = NORM_QC.out.pairs.map { f -> tuple(f.name.replace('.norm_qc_pairs.tsv', ''), f) }
                .join(wide_qc, failOnMismatch: true, failOnDuplicate: true)   // tuple(norm_run, pairs, wide table)
                .combine(tmm_counts)
                .map { nr, pairs, wide, bs, counts, libsizes -> tuple(nr, pairs, wide, counts, libsizes) }
            NORM_QC_PLOT(norm_plot_in, sc.norm_qc_plot)
            report_plots = report_plots.mix(NORM_QC_PLOT.out.svg)
        }
    }

    // ---- E. bigWigs, log2 ratio, TSS enrichment ----------------------------------
    def factors = Channel.empty()   // tuple(sample_id, multiplier, method)
    if (cfg.normRuns) {
        factors = sf_all.flatMap { nr, tsv ->
            tsv.splitCsv(header: true, sep: '\t').collect { r -> tuple(r.sample_id, r.bigwig_multiplier, r.method.toString().toLowerCase()) }
        }
        // IgG controls with no norm_run: multiplier 1, labelled 'unscaled' in the file name.
        def unscaled = cfg.activeIds.findAll { id -> !cfg.normRunOf[id] }
        if (unscaled) factors = factors.mix(Channel.fromList(unscaled.collect { id -> tuple(id, '1', 'unscaled') }))
    } else {
        factors = Channel.fromList(cfg.activeIds.collect { id -> tuple(id, '1', 'cpm') })
    }
    // final_bams holds active samples only, and every active sample is in exactly one
    // scale-factor table, so the keyed join must match one-to-one.
    def bw_in = final_bams.join(factors, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, bam, bai, multiplier, method)
    BIGWIG(bw_in, params.bigwig_bin_size, params.bigwig_units, isTrue(params.bigwig_extend_reads).toString(), params.effective_genome_size, sc.bigwig)
    def bw_lookup = fileLookup(BIGWIG.out.bw)   // id -> bigWig, every active sample

    def l2r_bw = Channel.empty()   // tuple(sample_id, log2ratio bigWig)
    if (isTrue(params.run_log2ratio) && cfg.controlOf) {
        def l2r_in = Channel.fromList(cfg.controlOf.collect { tid, cid -> tuple(tid, cid) })
            .combine(bw_lookup)
            .map { tid, cid, bk -> tuple(tid, pick(bk, tid, 'BIGWIG_LOG2RATIO'), cid, pick(bk, cid, 'BIGWIG_LOG2RATIO')) }
        BIGWIG_LOG2RATIO(l2r_in, params.bigwig_bin_size, sc.bigwig_log2ratio)
        l2r_bw = BIGWIG_LOG2RATIO.out.bw
    }

    TSS_ENRICHMENT(BIGWIG.out.bw, Channel.value(cfg.refs.gtf ?: cfg.refs.tss_bed), cfg.refs.gtf ? 'gtf' : 'bed', sc.qc_tss_enrichment)
    QC_TSS_TABLE(TSS_ENRICHMENT.out.tsv.collect(), samplesheet, sc.qc_tss_table)
    report_tables = report_tables.mix(QC_TSS_TABLE.out.table)
    report_plots = report_plots.mix(TSS_ENRICHMENT.out.svg)

    // ---- F. Peaks --------------------------------------------------------------
    def peak_files = Channel.empty()   // every peak BED, for QC_PEAK_STATS
    def frip_pairs = Channel.empty()   // tuple(sample_id, peaks_bed, peak_source)
    def csaw_by_name = Channel.empty()   // tuple(peakset, peaks bed)
    def class_by_name = Channel.empty()  // tuple(<id>.retained | <id>.lost | <id>.unclassified, bed)
    def annot_peaks = Channel.empty()   // tuple(label, bed)
    def pq = peakQcPlan(cfg)             // peak QC plan (run_peak_qc)
    def pq_sets = Channel.empty()        // tuple(peak set id, mark, bed): every peak set of the run

    if ((peakOn('run_csaw') || peakOn('run_csaw_cutoff_scan')) && cfg.peaksets) {
        def csaw_in = Channel.fromList(cfg.peaksets).combine(bam_lookup).map { ps, lk -> csawInput(ps, lk, 'CSAW') }
        if (peakOn('run_csaw')) {
            CSAW(csaw_in, sc.callpeaks_csaw)
            peak_files = peak_files.mix(CSAW.out.peaks)
            report_plots = report_plots.mix(CSAW.out.mds)
            csaw_by_name = CSAW.out.peaks.map { f -> tuple(f.name.replace('.peaks.bed', ''), f) }
            pq_sets = pq_sets.mix(csaw_by_name.map { id, f -> tuple(id.toString(), pq.sets[id.toString()].mark, f) })
            annot_peaks = annot_peaks.mix(csaw_by_name)
            // FRiP of each peak set in the samples of its target_group.
            def peakset_targets = cfg.peaksets.collectEntries { ps -> [(ps.peakset): ps.target_ids] }
            def csaw_frip = CSAW.out.peaks
                .map { f -> tuple(f.name.replace('.peaks.bed', ''), f) }
                .flatMap { pname, f -> peakset_targets[pname].collect { id -> tuple(id, f, "csaw_${pname}".toString()) } }
            frip_pairs = frip_pairs.mix(csaw_frip)
        }
        if (peakOn('run_csaw_cutoff_scan')) {
            def scan_in = Channel.fromList(cfg.peaksets).combine(bam_lookup).map { ps, lk -> csawInput(ps, lk, 'CSAW_CUTOFF_SCAN') }
            CSAW_CUTOFF_SCAN(scan_in, sc.csaw_cutoff_scan)
        }
    }

    // Target (non-control) members of every group that has any, in samplesheet order.
    def target_groups = cfg.groups
        .collect { g, ids -> [group: g, targets: ids.findAll { id -> !cfg.isControl.contains(id) }] }
        .findAll { gr -> gr.targets }

    if (cfg.runMacs3 && target_groups) {
        def macs_runs = target_groups.collect { gr ->
            [group: gr.group, targets: gr.targets, controls: gr.targets.collect { id -> cfg.controlOf[id] }.findAll { c -> c }.unique(), broad: cfg.macsBroad[gr.group].toString()]
        }
        def macs_in = Channel.fromList(macs_runs).combine(bam_lookup).map { gr, lk ->
            tuple(gr.group, bamsOf(lk, gr.targets, 'MACS3'), bamsOf(lk, gr.controls, 'MACS3 (control)'), gr.broad)
        }
        MACS3(macs_in, params.effective_genome_size, params.macs3_qvalue, params.macs3_broad_cutoff, sc.callpeaks_macs3)
        peak_files = peak_files.mix(MACS3.out.peaks.map { g, f -> f })
        annot_peaks = annot_peaks.mix(MACS3.out.peaks.map { g, f -> tuple("macs3_${g}".toString(), f) })
        pq_sets = pq_sets.mix(MACS3.out.peaks.map { g, f -> tuple("macs3_${g}".toString(), pq.sets["macs3_${g}".toString()].mark, f) })
        def macs_targets = target_groups.collectEntries { gr -> [(gr.group): gr.targets] }
        frip_pairs = frip_pairs.mix(MACS3.out.peaks.flatMap { g, f -> macs_targets[g].collect { id -> tuple(id, f, 'macs3') } })
    }

    if (cfg.runSeacr && cfg.controlOf) {
        SEACR_BEDGRAPH(final_bams, sc.seacr_bedgraph)
        def bg_lookup = fileLookup(SEACR_BEDGRAPH.out.bedgraph)
        def seacr_in = Channel.fromList(cfg.controlOf.collect { tid, cid -> tuple(tid, cid) })
            .combine(bg_lookup)
            .map { tid, cid, bk -> tuple(tid, pick(bk, tid, 'SEACR_CALL'), pick(bk, cid, 'SEACR_CALL (control)'), cfg.seacrMode[tid]) }
        SEACR_CALL(seacr_in, sc.seacr_call)
        def seacr_lookup = fileLookup(SEACR_CALL.out.bed)

        // Consensus per group: the called (non-control, with a control) members.
        def seacr_groups = target_groups
            .collect { gr -> [group: gr.group, called: gr.targets.findAll { id -> cfg.controlOf.containsKey(id) }] }
            .findAll { gr -> gr.called }
        def consensus_in = Channel.fromList(seacr_groups).combine(seacr_lookup).map { gr, sk ->
            tuple(gr.group, gr.called.collect { id -> pick(sk, id, 'SEACR_CONSENSUS') })
        }
        SEACR_CONSENSUS(consensus_in, params.seacr_min_replicates, params.seacr_reciprocal_overlap, sc.seacr_consensus)
        peak_files = peak_files.mix(SEACR_CALL.out.bed.map { sid, f -> f }, SEACR_CONSENSUS.out.bed.map { g, f -> f })
        annot_peaks = annot_peaks.mix(SEACR_CONSENSUS.out.bed.map { g, f -> tuple("seacr_${g}".toString(), f) })
        pq_sets = pq_sets.mix(SEACR_CONSENSUS.out.bed.map { g, f -> tuple("seacr_${g}".toString(), pq.sets["seacr_${g}".toString()].mark, f) })
        def seacr_members = seacr_groups.collectEntries { gr -> [(gr.group): gr.called] }
        frip_pairs = frip_pairs.mix(SEACR_CONSENSUS.out.bed.flatMap { g, f -> seacr_members[g].collect { id -> tuple(id, f, 'seacr') } })
    }

    // ---- F2. Region classification (run_classify) ------------------------------------------
    if (isTrue(params.run_classify) && cfg.classifications) {
        // tuple(classification, regions, ref ids, ref names, ref bams, ref bais, treat ids, treat names, treat bams, treat bais,
        //       norm_run, scale factors, lfc_retained, lfc_lost, fdr)
        def cls_in = Channel.fromList(cfg.classifications).combine(bam_lookup).combine(fileLookup(sf_all)).combine(fileLookup(csaw_by_name))
            .map { c, lk, sk, ck ->
                tuple(c.classification, bedLookupInput(c.regions, ck, 'CLASSIFY_REGIONS'),
                      c.ref_ids.join(','), namesOf(c.ref_ids, '.final.bam').join(','), bamsOf(lk, c.ref_ids, 'CLASSIFY_REGIONS'), baisOf(lk, c.ref_ids, 'CLASSIFY_REGIONS'),
                      c.treat_ids.join(','), namesOf(c.treat_ids, '.final.bam').join(','), bamsOf(lk, c.treat_ids, 'CLASSIFY_REGIONS'), baisOf(lk, c.treat_ids, 'CLASSIFY_REGIONS'),
                      c.norm_run, pick(sk, c.norm_run, 'CLASSIFY_REGIONS (scale factors)'), c.lfc_retained, c.lfc_lost, c.fdr)
            }
        CLASSIFY_REGIONS(cls_in, sc.classify_regions)
        class_by_name = CLASSIFY_REGIONS.out.beds.flatMap { id, r, l, u ->
            [tuple("${id}.retained".toString(), r), tuple("${id}.lost".toString(), l), tuple("${id}.unclassified".toString(), u)]
        }
        CLASSIFY_REGIONS_PLOT(CLASSIFY_REGIONS.out.per_region.join(CLASSIFY_REGIONS.out.summary, failOnMismatch: true, failOnDuplicate: true), sc.classify_regions_plot)
    }

    // ---- G. FRiP and peak statistics -------------------------------------------
    // One task per (sample, peak source). The sample's own flagstat comes from a map
    // keyed by sample id (a join cannot be used: a sample can have several sources).
    def any_peaks = (peakOn('run_csaw') && cfg.peaksets) || (cfg.runMacs3 && target_groups) || (cfg.runSeacr && cfg.controlOf)
    if (any_peaks) {
        def frip_in = frip_pairs.combine(bam_lookup).combine(flagstat_lookup).map { sid, peaks, source, lk, fk ->
            tuple(sid, pick(lk, sid, 'FRIP')[0], pick(lk, sid, 'FRIP')[1], peaks, source, pick(fk, sid, 'FRIP (flagstat)'))
        }
        FRIP(frip_in, sc.qc_frip)
        QC_FRIP_TABLE(FRIP.out.tsv.collect(), samplesheet, sc.qc_frip_table)
        report_tables = report_tables.mix(QC_FRIP_TABLE.out.table)

        QC_PEAK_STATS(peak_files.collect(), sc.qc_peak_stats)
        report_tables = report_tables.mix(QC_PEAK_STATS.out.table)
    }

    // csaw peak sets and classification regions by name (one lookup for every region reference)
    def base_lk = fileLookup(csaw_by_name.mix(class_by_name))
    def intersect_by_name = Channel.empty()   // tuple(<intersect name>, intersect bed): intersect names are regions too

    // ---- I4. Intersections (before the heatmaps, tracks and annotation, which may name an intersect) -----
    if (isTrue(params.run_intersect) && cfg.intersects) {
        def inter_in = Channel.fromList(cfg.intersects).combine(base_lk)
            .map { r, nk -> tuple(r.name, r.label_a, r.label_b, r.min_overlap ?: '', r.reciprocal, r.output, bedLookupInput(r.bed_a, nk, 'INTERSECT'), bedLookupInput(r.bed_b, nk, 'INTERSECT')) }
        INTERSECT(inter_in, sc.intersect_beds)
        intersect_by_name = INTERSECT.out.bed.map { f -> tuple(f.name.replace('.intersect.bed', ''), f) }
    }

    // one lookup for every region reference: csaw peak sets, classification regions and intersects
    def named_lk = fileLookup(csaw_by_name.mix(class_by_name, intersect_by_name))

    // ---- I. Heatmaps -----------------------------------------------------------------
    if (isTrue(params.run_heatmaps) && cfg.heatmaps) {
        // tuple(id, hm, region files, bigwig list); regions = one file per region of the row
        def hm_in = Channel.fromList(cfg.heatmaps).combine(bw_lookup).combine(named_lk)
            .map { hm, bk, nk -> tuple(hm.heatmap, hm, hm.regions.collect { r -> bedLookupInput(r, nk, 'HEATMAPS') }, hm.samples.collect { id -> pick(bk, id, 'HEATMAPS') }) }
        def hm_dt = hm_in.filter { id, hm, regions, bws -> hm.engine == 'deeptools' }
        def hm_eh = hm_in.filter { id, hm, regions, bws -> hm.engine == 'enrichedheatmap' }

        // deepTools: one matrix per row (unclustered rows are plotted from it)
        COMPUTEMATRIX(hm_dt.map { id, hm, regions, bws -> tuple(id, id, regions.collect { f -> f.name }.join(','), regions, hm.samples.join(','), bws.collect { f -> f.name }.join(','), bws, hm.window_bp, hm.bin_bp) }, sc.heatmap_computematrix)
        def dt_mat = COMPUTEMATRIX.out.matrix.join(hm_dt, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, matrix, hm, regions, bws)

        // EnrichedHeatmap: one matrix per row, plus the cluster matrix of the clustered rows (one region only)
        BUILD_MATRIX(hm_eh.map { id, hm, regions, bws -> tuple(id, regions[0], hm.samples.join(','), bws.collect { f -> f.name }.join(','), bws, hm.window_bp, hm.bin_bp, hm.cluster_k ? hm.cluster_samples.join(',') : '') }, sc.heatmap_build_matrix)
        def eh_mat = BUILD_MATRIX.out.mats.join(hm_eh, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, rds, targets, hm, regions, bws)

        // Clustering of both engines (only rows with cluster_k, which take one region). The cluster matrix exists only for
        // clustered EnrichedHeatmap rows, so those rows alone are joined with it.
        def cl_dt = dt_mat.filter { id, mat, hm, regions, bws -> hm.cluster_k }
            .map { id, mat, hm, regions, bws -> tuple(id, mat, regions[0], hm.cluster_samples.join(','), hm.cluster_k) }
        def cl_eh = BUILD_MATRIX.out.cluster_matrix
            .join(hm_eh.filter { id, hm, regions, bws -> hm.cluster_k }, failOnMismatch: true, failOnDuplicate: true)
            .map { id, cm, hm, regions, bws -> tuple(id, cm, regions[0], '', hm.cluster_k) }
        HEATMAP_CLUSTER(cl_dt.mix(cl_eh), params.heatmap_seed, sc.heatmap_cluster)

        // deepTools clustered: the matrix is recomputed with the cluster BEDs as separate region groups
        // (never plotHeatmap --kmeans). The plain join keeps only deepTools ids: EnrichedHeatmap ids
        // have no deepTools row, so they drop out here on purpose.
        def cl_beds = HEATMAP_CLUSTER.out.beds.map { id, beds -> tuple(id, kSorted(beds instanceof List ? beds : [beds])) }
        def hm_dt_cl = hm_dt.filter { id, hm, regions, bws -> hm.cluster_k }
        def cl_dt_in = cl_beds.join(hm_dt_cl)   // tuple(id, beds, hm, regions, bws)
        COMPUTEMATRIX_CLUSTERED(cl_dt_in.map { id, beds, hm, regions, bws -> tuple(id, "${id}.clustered".toString(), beds.collect { b -> b.name }.join(','), beds, hm.samples.join(','), bws.collect { f -> f.name }.join(','), bws, hm.window_bp, hm.bin_bp) }, sc.heatmap_computematrix)
        def cl_labels = cl_dt_in.map { id, beds, hm, regions, bws -> tuple(id, regionLabels(hm.cluster_labels, beds.size()).join(',')) }
        def cl_label_lk = fileLookup(cl_labels)   // id -> region labels, one value

        // deepTools plots: unclustered matrices from COMPUTEMATRIX, clustered ones from COMPUTEMATRIX_CLUSTERED
        def plain_mat = dt_mat.filter { id, mat, hm, regions, bws -> !hm.cluster_k }.map { id, mat, hm, regions, bws -> tuple(id, mat, hm) }
        def clustered_mat = COMPUTEMATRIX_CLUSTERED.out.matrix.join(hm_dt_cl.map { id, hm, regions, bws -> tuple(id, hm) }, failOnMismatch: true, failOnDuplicate: true)   // tuple(id, matrix, hm)

        // Colour limits of 'auto' z_max / y_max: HEATMAP_LIMITS on the plotted matrix, one limit per sample.
        def hm_limits_in = plain_mat.mix(clustered_mat).filter { id, mat, hm -> hm.auto }
            .map { id, mat, hm -> tuple(id, mat, hm.samples.join(','), hm.antibodies.join(',')) }
        HEATMAP_LIMITS(hm_limits_in, sc.heatmap_limits)
        def lim_lk = fileLookup(HEATMAP_LIMITS.out.limits.map { id, f -> tuple(id, f.splitCsv(sep: '\t', header: true).collect { r -> r.limit }) })   // id -> limits, sample order

        def plot_in = plain_mat.mix(clustered_mat).combine(cl_label_lk).combine(lim_lk)
            .map { id, mat, hm, cl, lim ->
                def n = hm.samples.size()
                def lims = lim.containsKey(id) ? lim[id] : []
                tuple(id, mat, hm.labels.join(','), hm.colors.join(','), hm.z_min.join(','),
                      resolveAuto(hm.z_max, n, lims, "heatmap ${id} z_max").join(','),
                      resolveAuto(hm.y_max, n, lims, "heatmap ${id} y_max").join(','),
                      cl.containsKey(id) ? cl[id] : hm.region_labels.join(','), hm.sort_idx)
            }
        PLOTHEATMAP(plot_in, (params.heatmap_y_label ?: "normalised ${params.bigwig_units}").toString(), sc.heatmap_plot)

        // EnrichedHeatmap plots: cluster assignments looked up by id (no remainder joins)
        def asg_lk = fileLookup(HEATMAP_CLUSTER.out.assignments)
        PLOT_ENRICHED(eh_mat.combine(asg_lk).map { id, rds, tg, hm, regions, bws, lk -> tuple(id, rds, tg, hm.samples.join(','), hm.labels.join(','), hm.colors.join(','), hm.legend_titles.join(','), hm.z_min.join(','), hm.z_max.join(','), hm.y_max.join(','), hm.window_bp, hm.cluster_labels.join(','), hm.sort_using, lk.containsKey(id) ? lk[id] : []) }, sc.heatmap_plot_enriched)
    }

    // ---- I1b. Peak QC regions (run_peak_qc) --------------------------------------------------------------
    // Automatic QC loci chosen from the peak sets of this run (bin/peak_qc_regions.py). With run_tracks they are drawn as the
    // track set peak_qc by the tracks chain below. Skipped silently when the run makes no peak set.
    def peak_qc_on = peakOn('run_peak_qc') && pq.sets
    if (peak_qc_on && !params.gtf) {
        log.warn "run_peak_qc: gtf is not set in config.yaml: peak QC regions are skipped"
        peak_qc_on = false
    }
    if (peak_qc_on) {
        def pq_sets_tsv = pq_sets.map { id, mark, f -> "${id}\t${mark}\t${f.name}".toString() }
            .collectFile(name: 'peaksets.tsv', newLine: true, seed: 'id\tmark\tbed', sort: true)
        def pq_cls_tsv = Channel.value([])
        def pq_cls_tables = Channel.empty()
        if (pq.classes) {
            def pq_pr = CLASSIFY_REGIONS.out.per_region.filter { id, f -> pq.classes.containsKey(id.toString()) }
            pq_cls_tsv = pq_pr.map { id, f ->
                    def c = pq.classes[id.toString()]
                    "${id}\t${c.mark}\t${c.base}\t${f.name}\t${c.contrast}\t${c.ref_group}".toString()
                }
                .collectFile(name: 'classes.tsv', newLine: true, seed: 'id\tmark\tbase_peakset\ttable\tcontrast\tref_group', sort: true)
            pq_cls_tables = pq_pr.map { id, f -> f }
        }
        def pq_bw_tsv = bw_lookup.flatMap { bk ->
                pq.bw_ids.collect { id -> "${id}\t${pq.mark_of[id]}\t${pick(bk, id, 'PEAK_QC_REGIONS').name}\t${pq.group_of[id]}".toString() }
            }
            .collectFile(name: 'bigwigs.tsv', newLine: true, seed: 'sample\tmark\tbigwig\tcondition', sort: true)
        def pq_bw_files = bw_lookup.flatMap { bk -> pq.bw_ids.collect { id -> pick(bk, id, 'PEAK_QC_REGIONS') } }
        def pq_files = pq_sets.map { id, mark, f -> f }.mix(pq_cls_tables, pq_bw_files).toList()
        def pq_excl = csvList(params.peak_qc_exclude_bed).collect { b -> file(hostPath(b)) }
        def pq_shape = textOf(params.peak_qc_shape_marks) ?: (pq.marks.size() > 3 ? pq.marks.take(3).join(',') : '')
        // IgG cap of the nopeak_highsignal tiles: only passed when an IgG bigWig is among the peak-QC bigWigs (the script stops otherwise)
        def pq_igg_cap = pq.bw_ids.any { id -> pq.mark_of[id].equalsIgnoreCase('IgG') } ? textOf(params.peak_qc_nopeak_igg_max_fold) : ''
        PEAK_QC_REGIONS(pq_sets_tsv, pq_cls_tsv, pq_bw_tsv, pq_files, Channel.value(cfg.refs.gtf), Channel.value(cfg.refs.chrom_sizes),
                        Channel.value(cfg.refs.blacklist ?: []), Channel.value(pq_excl), pq_excl.collect { b -> b.name }.join(','), pq_shape, pq_igg_cap, sc.peak_qc_regions)
    }

    // ---- I2. Tracks -------------------------------------------------------------------
    if (isTrue(params.run_tracks) && (cfg.tracks || peak_qc_on)) {
        // tuple(track_set, track, bigwigs, bed files, log2 sample ids, log2 bigWigs)
        def trk = Channel.fromList(cfg.tracks).combine(bw_lookup).combine(fileLookup(l2r_bw)).combine(named_lk)
            .map { t, bk, lk, nk ->
                def bws = t.samples.collect { id -> pick(bk, id, 'TRACKS') }
                def beds = t.beds.collect { b -> bedLookupInput(b, nk, 'TRACKS') }
                def l2ids = t.log2ratio ? t.samples.findAll { id -> lk.containsKey(id) } : []
                def l2 = l2ids.collect { id -> lk[id] }
                tuple(t.track_set, t, bws, beds, l2ids, l2)
            }
        if (peak_qc_on) {
            // Automatic track set peak_qc: bigWigs of the peak sets' target groups (+ IgG), y max shared per antibody
            // (tracks_max_value_default), BED tracks = every peak set + the retained / lost regions of the classifications.
            // Peak-set BED tracks: 'csaw <peakset>' black, 'MACS3 <group>' blue, 'SEACR <group>' vermillion (macs3_ / seacr_ prefix stripped from the id).
            def pq_bed_ch = pq_sets.map { id, mark, f ->
                def sid = id.toString()
                sid.startsWith('macs3_') ? tuple(sid, f, '#0072B2', "MACS3 ${sid.substring(6)}".toString())
                    : sid.startsWith('seacr_') ? tuple(sid, f, '#D55E00', "SEACR ${sid.substring(6)}".toString())
                    : tuple(sid, f, '#000000', "csaw ${sid}".toString())
            }
            if (pq.classes) {
                pq_bed_ch = pq_bed_ch.mix(class_by_name
                    .filter { n, f -> !n.endsWith('.unclassified') && pq.classes.containsKey(n.substring(0, n.lastIndexOf('.'))) }
                    .map { n, f -> tuple(n.toString(), f, n.endsWith('.retained') ? '#1f78b4' : '#e31a1c', n.replace('.', ' ').toString()) })
            }
            def pq_beds = pq_bed_ch.toSortedList { a, b -> a[0] <=> b[0] }.map { rows -> [rows: rows] }
            def pq_trk = PEAK_QC_REGIONS.out.regions
                .filter { f -> f.readLines().any { l -> l.trim() } }
                .combine(bw_lookup).combine(pq_beds)
                .map { regions_f, bk, bd ->
                    def loci = regions_f.readLines().findAll { l -> l.trim() }.collect { l -> l.split('\t') }.collect { r -> "${r[0]}=${r[1]}".toString() }
                    def t = [track_set: 'peak_qc', loci: loci, flanks: ['0'], samples: pq.track_ids, labels: pq.track_ids, colors: pq.track_colors,
                             max_tokens: pq.track_max, bed_colors: bd.rows.collect { b -> b[2] }, bed_labels: bd.rows.collect { b -> b[3] }, scalebar: 'auto', log2ratio: false,
                             genes_fontsize: 7, genes_height: 1.5, scalebar_fontsize: 12]
                    tuple('peak_qc', t, pq.track_ids.collect { id -> pick(bk, id, 'TRACKS (peak_qc)') }, bd.rows.collect { b -> b[1] }, [], [])
                }
            trk = trk.mix(pq_trk)
        }
        TRACKS_INI(trk.map { ts, t, bws, beds, l2ids, l2 ->
                tuple(ts, t.loci.join(','), t.flanks.join(','), t.samples.join(','), bws.collect { f -> f.name }.join(','), t.labels.join(','), t.colors.join(','), t.max_tokens.join(','),
                      (0..<beds.size()).collect { i -> stageBedName(i, beds[i]) }.join(','), t.bed_colors.join(','), t.bed_labels.join(','), t.scalebar, l2.collect { f -> f.name }.join(','), l2ids.join(','),
                      [t.genes_fontsize != null ? "--genes-fontsize=${t.genes_fontsize}" : null, t.genes_height != null ? "--genes-height=${t.genes_height}" : null, t.scalebar_fontsize != null ? "--scalebar-fontsize=${t.scalebar_fontsize}" : null, t.ini_args ? t.ini_args : null].findAll { x -> x }.join(' ').toString())
            }, Channel.value(cfg.refs.gtf), sc.make_tracks_ini)
        // one pyGenomeTracks figure per region of each track set; the limits of group:NAME tracks are resolved per locus first
        def ini_regions = TRACKS_INI.out.ini.flatMap { ts, ini, tsv -> tsv.splitCsv(sep: '\t').collect { r -> tuple(ts, r[0], r[1], ini) } }
        def trk_files = trk.map { ts, t, bws, beds, l2ids, l2 -> tuple(ts, bws, beds, l2) }
        TRACKS_LIMITS(ini_regions.combine(trk_files.map { ts, bws, beds, l2 -> tuple(ts, bws) }, by: 0), params.tracks_limit_stat, sc.tracks_limits)
        PYGENOMETRACKS(TRACKS_LIMITS.out.ini.combine(trk_files, by: 0), Channel.value(cfg.refs.gtf), params.tracks_width, params.tracks_dpi, sc.tracks_pygenometracks)
    }

    // ---- I3. Peak annotation ----------------------------------------------------------
    if (peakOn('run_annotation') && cfg.annotation) {
        annot_peaks = annot_peaks.mix(Channel.fromList(cfg.annotationExtra).combine(named_lk).map { x, nk -> tuple(x.label, bedLookupInput(x.bm, nk, 'PEAK_ANNOTATION')) })
        def annot_in = annot_peaks.toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows -> tuple(rows.collect { r -> r[0] }.join(','), rows.collect { r -> r[1].name }.join(','), rows.collect { r -> r[1] }) }
        PEAK_ANNOTATION(annot_in, Channel.value(cfg.refs.gtf), params.annotation_colors ?: '', isTrue(params.annotation_kegg).toString(), params.annotation_tss_region, sc.peak_annotation)
    }

    // ---- H. Report ----------------------------------------------------------------
    if (isTrue(params.run_report)) {
        REPORT(report_tables.collect(), report_plots.collect(), params.project_name, sc.build_report)
    }
}
