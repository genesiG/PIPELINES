#!/usr/bin/env nextflow
/*
 * CRISPR_SCREEN — reproducible CRISPR screen pipeline
 * ======================================================
 * Demultiplex (or take already-demultiplexed FASTQs) -> FastQC/MultiQC ->
 * mageck count per library -> guide-representation QC -> (optional) mageck
 * test, fanned out per contrast.
 *
 * Docker version (Windows / Mac): run from your own project directory with
 *   run_crisprscreen
 * which starts the pipeline's Docker image and runs, inside it,
 *   nextflow run /pipeline -params-file config.yaml -resume
 * (docker/run_in_container.sh). Every step runs in that same container.
 *
 * See README.md for the samplesheet.csv / library.csv / contrasts.csv schemas.
 */

nextflow.enable.dsl = 2

include { DEMULTIPLEX_SE; DEMULTIPLEX_PE } from './modules/local/demultiplex'
include { PLOT_BARCODE_POSITION }           from './modules/local/plot_barcode_position'
include { FASTQC }                          from './modules/local/fastqc'
include { MULTIQC }                          from './modules/local/multiqc'
include { CUTADAPT_STATS }                    from './modules/local/cutadapt_stats'
include { BUILD_MAGECK_LIBRARY }               from './modules/local/build_mageck_library'
include { GUIDE_ANCHOR_DISCOVERY }               from './modules/local/guide_anchor_discovery'
include { GUIDE_ANCHOR_TRIM }                     from './modules/local/guide_anchor_trim'
include { MAGECK_COUNT }                        from './modules/local/mageck_count'
include { LIBRARY_BREAKDOWN }                   from './modules/local/library_breakdown'
include { GUIDE_REPRESENTATION_WRANGLE; GUIDE_REPRESENTATION_PLOT } from './modules/local/guide_representation'
include { GUIDE_REPRESENTATION_REPORT }                              from './modules/local/guide_representation_report'
include { GUIDE_REPRESENTATION_COMPARE_REPORT }                       from './modules/local/guide_representation_compare_report'
include { VALIDATE_CONTRASTS }                    from './modules/local/validate_contrasts'
include { MAGECK_TEST }                             from './modules/local/mageck_test'
include { SOFTWARE_VERSIONS }                        from './modules/local/software_versions'
include { METHODS }                                   from './modules/local/methods'

// -----------------------------------------------------------------------
// Parameter validation + samplesheet parsing, done once at run start.
// Nextflow's strict-syntax parser forbids bare top-level statements
// (if/for/etc.) outside a process, workflow or function — so all of this
// setup logic lives in one function, called from the top of workflow{}.
// -----------------------------------------------------------------------
// Every step runs a script from this repo's bin/. A project replaces one by
// setting a config.yaml key named after the script (its file name without
// the extension) to the path of its own copy, e.g.
//     mageck_count: utils/mageck_count.sh
// Relative paths resolve against the project directory. Returns
// [script name: file] for every script in bin/, replacements applied.
def stepScripts() {
    def defaults = defaultScripts()

    // Keys from before scripts were replaced by name.
    def renamed = [plot_theme: 'my_theme', barcode_position_plot_script: 'plot_barcode_position',
                   barcode_position_plot_args: 'plot_barcode_position_args']
    ['barcode_position_qc', 'cutadapt_stats', 'build_mageck_library', 'validate_contrasts',
     'guide_representation_wrangle', 'guide_representation_plot'].each { n -> renamed["${n}_script".toString()] = n }
    def old = renamed.findAll { k, v -> params.containsKey(k) && params[k] }.collect { k, v -> "${k} -> ${v}" }
    if (old) fail("these config.yaml keys were renamed (a script is now replaced by a key named after it):\n  " + old.join('\n  '))

    // A script path under a key that names no pipeline script is most likely a typo.
    def stray = params.findAll { k, v -> v instanceof CharSequence && v.toString() =~ /(?i)\.(r|py|sh)$/ && !defaults.containsKey(k) }
    if (stray) {
        fail("config.yaml: no pipeline script is named ${stray.keySet().join(', ')}. To replace a step, name the key " +
             "after the script's file name without its extension. Scripts: ${defaults.keySet().sort().join(', ')}")
    }

    def scripts = [:]
    def record = readDefaultsRecord()
    defaults.each { name, f ->
        def replacement = params.containsKey(name) ? params[name] : null
        if (replacement) {
            def r = file(replacement.toString())
            if (!r.exists()) fail("config.yaml: ${name}: ${replacement} not found.")
            def dflt = f
            def rel = projectPath(r)
            def recorded = record[rel]
            def current = md5NoCR(r)
            def current_default = md5NoCR(dflt)
            // The copy was made from an older default: refresh it if unedited, otherwise say so.
            if (recorded && recorded != current_default) {
                if (current == recorded) {
                    try {
                        java.nio.file.Files.write(r, dflt.bytes)
                        updateDefaultsRecord(rel, current_default)
                        log.info "${rel} was an unedited copy of the pipeline's default and has been updated to the current version"
                    } catch (Exception e) {
                        log.warn "${rel}: could not update it to the pipeline's current ${dflt.name} (${e.message}); keeping the project's copy"
                    }
                } else {
                    log.warn "${rel}: the pipeline's default ${dflt.name} changed since this copy was made; this project's edited copy is used. Compare it with ${projectDir}/bin/${dflt.name} (the copy was made from an older version)"
                }
            }
            def note = md5NoCR(r) == current_default ? " (an unchanged copy of the pipeline's default)" : ''
            log.info "Using this project's ${name}: ${replacement}${note}"
            scripts[name] = r
        } else {
            scripts[name] = f
        }
    }
    return scripts
}

// Config keys that were removed: say so instead of silently ignoring them.
def warnRemovedKeys() {
    def removed = [target_mismatches: 'it only named output files',
                   mageck_use_r2    : 'it was never used']
    removed.each { key, why ->
        if (params.containsKey(key)) log.warn "${key} is no longer used (${why}) — remove it from config.yaml"
    }
}

def fail(String message) {
    log.error "ERROR: ${message}"
    System.exit(1)
}

// Container image of every process: SOFTWARE_VERSIONS probes each distinct
// image once, and METHODS.md says which image ran which step. Docker
// version: every process runs in the one image this Nextflow runs in
// (CRISPR_IMAGE, set by run_crisprscreen). Returns [process name: image].
def processContainers() {
    def image = System.getenv('CRISPR_IMAGE') ?: 'crisprscreen'
    def images = [:]
    file("${projectDir}/modules/local").listFiles().findAll { it.name.endsWith('.nf') }.sort { it.name }.each { f ->
        f.readLines().each { line ->
            def p = line =~ '^\\s*process\\s+(\\w+)\\s*\\{'
            if (p.find()) images[p.group(1)] = image
        }
    }
    return images
}

def md5(f) {
    return java.security.MessageDigest.getInstance('MD5').digest(f.bytes).encodeHex().toString()
}

// md5 of the file's bytes without CR (13) bytes: a copy that git checked out with CRLF
// line endings (core.autocrlf on Windows) compares equal to the pipeline's LF original.
def md5NoCR(f) {
    def kept = new java.io.ByteArrayOutputStream()
    f.bytes.each { b -> if (b != 13) kept.write(b as int) }
    return java.security.MessageDigest.getInstance('MD5').digest(kept.toByteArray()).encodeHex().toString()
}

// Path relative to the project (launch) directory when it lies inside it.
def projectPath(f) {
    def p = f.toString()
    def root = "${workflow.launchDir}/".toString()
    return p.startsWith(root) ? p.substring(root.length()) : p
}

// The input path as the user wrote it on this computer. The Docker version's params
// hold container paths and run_crisprscreen passes the host path in the environment;
// on HPC the variable is unset, so the param (already the host path) is used.
def hostInput(String envName, def param) {
    def h = System.getenv(envName)
    return h ? h : (param ? param.toString() : null)
}

// [script name: file] for every script in bin/ (the pipeline's defaults).
def defaultScripts() {
    def defaults = [:]
    file("${workflow.projectDir}/bin").listFiles().findAll { it.isFile() }.each { f -> defaults[f.baseName] = f }
    return defaults
}

// utils/.pipeline_defaults.md5 (written by `init`): md5 of each default script copied into the project, as copied. [project-relative path: md5]
def defaultsRecordFile() { return file("${workflow.launchDir}/utils/.pipeline_defaults.md5") }
def readDefaultsRecord() {
    def out = [:]
    def rec = defaultsRecordFile()
    if (!rec.exists()) return out
    rec.readLines().each { line ->
        def m = line =~ /^([0-9a-fA-F]{32})\s+(\S.*?)\s*$/
        if (m.find()) out[m.group(2)] = m.group(1).toLowerCase()
    }
    return out
}
// Replaces the record line for path with newMd5; every other line is kept as it is.
def updateDefaultsRecord(String path, String newMd5) {
    def rec = defaultsRecordFile()
    if (!rec.exists()) return
    def newLine = "${newMd5}  ${path}".toString()
    def lines = rec.readLines().collect { line ->
        def m = line =~ /^([0-9a-fA-F]{32})\s+(\S.*?)\s*$/
        return (m.find() && m.group(2) == path) ? newLine : line
    }
    if (!lines.contains(newLine)) lines << newLine
    rec.text = lines.join('\n') + '\n'
}

// Trimmed stdout of an external command, or null if it fails or prints nothing.
def commandOutput(List cmd) {
    try {
        def proc = cmd.execute()
        def out = new StringBuilder()
        def err = new StringBuilder()
        proc.waitForProcessOutput(out, err)
        def text = out.toString().trim()
        return proc.exitValue() == 0 && text ? text : null
    } catch (Exception e) {
        return null
    }
}

// git in the pipeline repository (safe.directory: labmates run a repo they don't own).
def pipelineGit(List args) {
    return commandOutput(['git', '-c', 'safe.directory=*', '-C', projectDir.toString()] + args)
}

// [path, md5] of a small input table (samplesheet, library, contrasts), or null if unset.
def inputFile(p) {
    return p ? [path: projectPath(file(p.toString())), md5: md5(file(p.toString()))] : null
}

// This run's resolved configuration, for METHODS.md (bin/write_methods.py):
// what was run, on what, with which scripts and arguments. Facts the steps
// themselves record (tool command lines, versions, anchors) are read from
// their outputs instead.
def runInfo(cfg, boolean autoAnchor) {
    // Docker version: the pipeline folder is mounted read-only and the image
    // has no view of its .git, so run_crisprscreen reads the commit (and
    // whether this folder has uncommitted changes) on the host and passes it.
    def envCommit = System.getenv('CRISPR_PIPELINE_COMMIT')
    def commit = envCommit ?: pipelineGit(['rev-parse', '--short', 'HEAD'])
    def modified = envCommit ? System.getenv('CRISPR_PIPELINE_MODIFIED') == 'true'
                             : (commit ? pipelineGit(['status', '--porcelain', '--untracked-files=no']) != null : null)
    def scripts = [:]
    def defaults = defaultScripts()
    cfg.scripts.each { name, f ->
        def dflt = defaults[name]
        def replaced = !f.toString().startsWith("${projectDir}/bin/".toString())
        // A copy that is byte-identical to the pipeline's default is reported as the default.
        if (replaced && dflt != null && md5NoCR(f) == md5NoCR(dflt)) replaced = false
        scripts[name] = [path: replaced ? projectPath(f) : "bin/${(dflt != null ? dflt : f).name}".toString(), replaced: replaced, md5: md5(f)]
    }
    def extraArgs = [:]
    cfg.scripts.keySet().each { name ->
        def key = "${name}_args".toString()
        if (params.containsKey(key) && params[key]) extraArgs[name] = params[key].toString()
    }
    def simpleParams = [:]
    params.each { k, v ->
        if (v == null || v instanceof Boolean || v instanceof Number) simpleParams[k.toString()] = v
        else if (v instanceof CharSequence) simpleParams[k.toString()] = v.toString()
    }
    def paramsFile = workflow.commandLine =~ '-params-file\\s+(\\S+)'
    // Docker version: Nextflow runs inside the container (no engine of its
    // own), so the host's Docker version, the image and the host OS come from
    // run_crisprscreen via the environment.
    def engine = System.getenv('CRISPR_CONTAINER_ENGINE') ?: workflow.containerEngine
    def engineVersion = System.getenv('CRISPR_CONTAINER_ENGINE_VERSION') ?: (engine ? commandOutput([engine, '--version']) : null)
    return [
        pipeline : [name: workflow.manifest.name, version: workflow.manifest.version,
                    path: System.getenv('CRISPR_PIPELINE_HOST_DIR') ?: projectDir.toString(),
                    commit: commit, modified: modified],
        nextflow : [version: nextflow.version.toString(), run_name: workflow.runName, session_id: workflow.sessionId.toString(),
                    start: workflow.start.toString(), profile: workflow.profile, stub_run: workflow.stubRun,
                    launch_dir: workflow.launchDir.toString(), config_file: paramsFile.find() ? paramsFile.group(1) : null,
                    container_engine: engine, container_engine_version: engineVersion,
                    image: System.getenv('CRISPR_IMAGE'), image_id: System.getenv('CRISPR_IMAGE_ID'),
                    host: System.getenv('CRISPR_HOST_OS'),
                    cpus: System.getenv('CRISPR_MAX_CPUS'), memory_gb: System.getenv('CRISPR_MAX_MEMORY_GB')],
        params   : simpleParams,
        input    : [mode: cfg.alreadyDemuxed ? 'fastq_dir' : 'demultiplex', paired_end: params.paired_end as boolean,
                    fastq_r1: cfg.alreadyDemuxed ? null : hostInput('CRISPR_FASTQ_R1_HOST', params.fastq_r1),
                    fastq_r2: cfg.alreadyDemuxed ? null : hostInput('CRISPR_FASTQ_R2_HOST', params.fastq_r2),
                    fastq_dir: hostInput('CRISPR_FASTQ_DIR_HOST', params.fastq_dir)],
        samplesheet : inputFile(params.samplesheet),
        samples  : cfg.sampleNames,
        barcodes : cfg.barcodes,
        libraries: cfg.libraryGroups.collect { g -> [name: g.name, path: projectPath(g.file), md5: md5(g.file), samples: g.samples] },
        counting_mode : cfg.allVsAll ? 'all_vs_all' : (cfg.libraryGroups.size() > 1 ? 'assigned' : 'single'),
        guide_position: autoAnchor ? [mode: 'anchor'] : [mode: 'fixed', trim5: params.trim5.toString()],
        contrasts: inputFile(params.contrasts),
        scripts  : scripts,
        extra_args: extraArgs,
        containers: cfg.images
    ]
}

def prepareRun() {
    warnRemovedKeys()
    if (params.help) {
        log.info """
        CRISPR_SCREEN pipeline
        -----------------------
        Usage (from your project folder, which holds config.yaml):
          run_crisprscreen

        Required params (set in config.yaml):
          project_name    Prefix used for every output file
          samplesheet     CSV: sample_name (+ barcode_id,barcode when demultiplexing; + optional library)
          library         library CSV(s) (glob ok): gene,gRNA_ID,sequence
          and one of:
          fastq_dir       directory of already-demultiplexed <sample_name>.fastq.gz files —
                          demultiplexing is skipped (fastq_r1/fastq_r2 are then ignored)
          fastq_r1        multiplexed R1 (or single-end) fastq.gz — the pipeline demultiplexes it
                          (+ fastq_r2 if paired_end: true)

        Optional:
          contrasts       CSV: treatment,control — omit to stop after counting/QC
          outdir          Default: ./results
        """
        System.exit(0)
    }

    // ---- Required params / input mode --------------------------------------
    def missing = []
    if (!params.samplesheet) missing << 'samplesheet'
    if (!params.library)     missing << 'library'
    if (missing) fail("missing required param(s): ${missing.join(', ')} (set in your config.yaml).")

    // fastq_dir set = already demultiplexed; fastq_r1/fastq_r2 are then ignored.
    def alreadyDemuxed = params.fastq_dir as boolean
    if (!params.fastq_r1 && !params.fastq_dir) {
        fail("config.yaml needs fastq_r1 (a multiplexed FASTQ to demultiplex) or fastq_dir " +
             "(a directory of already-demultiplexed per-sample FASTQs).")
    }
    if (!alreadyDemuxed && params.paired_end && !params.fastq_r2) {
        fail("paired_end: true needs fastq_r2 (the multiplexed R2 FASTQ).")
    }

    // ---- Samplesheet --------------------------------------------------------
    def sheet = file(params.samplesheet)
    if (!sheet.exists()) fail("samplesheet not found: ${params.samplesheet}")
    // Normalize: strip a UTF-8 BOM and whitespace/CR from headers and values;
    // drop blank rows (e.g. a trailing ",,").
    def rows = sheet.splitCsv(header: true)
        .collect { r -> r.collectEntries { k, v -> [(k.toString().replace('﻿', '').trim()), (v ?: '').toString().trim()] } }
        .findAll { it.sample_name }
    if (!rows) fail("${params.samplesheet} has no rows with a sample_name.")
    def columns = rows[0].keySet()
    if (!('sample_name' in columns)) fail("${params.samplesheet} needs a sample_name column.")

    def sampleNames = rows.collect { it.sample_name }
    def dups = sampleNames.findAll { n -> sampleNames.count(n) > 1 }.unique()
    if (dups) fail("${params.samplesheet} has duplicate sample_name values: ${dups.join(', ')}")

    // Barcodes: only needed when the pipeline demultiplexes.
    def adapterArgs = null
    def barcodeLen  = null
    if (!alreadyDemuxed) {
        if (!('barcode' in columns)) {
            fail("${params.samplesheet} needs a barcode column to demultiplex fastq_r1. " +
                 "If your samples are already demultiplexed, set fastq_dir instead of fastq_r1.")
        }
        def noBarcode = rows.findAll { !it.barcode }.collect { it.sample_name }
        if (noBarcode) fail("these samples have no barcode in ${params.samplesheet}: ${noBarcode.join(', ')}")
        def barcodes = rows.collect { it.barcode }
        def dupBarcodes = barcodes.findAll { b -> barcodes.count(b) > 1 }.unique()
        if (dupBarcodes) {
            def collided = rows.findAll { it.barcode in dupBarcodes }
                .collect { "${it.sample_name} (${it.barcode})" }.join(', ')
            fail("${params.samplesheet} has samples sharing the same barcode — cutadapt would " +
                 "silently assign all matching reads to whichever sample is listed first and " +
                 "leave the rest empty: ${collided}")
        }
        adapterArgs = rows.collect { "-g ${it.sample_name}=${it.barcode}" }.join(' ')
        barcodeLen  = rows.collect { it.barcode.length() }.max()
    }

    // ---- Per-sample FASTQs ---------------------------------------------------
    // r1Name: sample -> the FASTQ file name mageck count reads for it.
    def r1Name        = [:]
    def sampleFastqs  = []   // fastq_dir mode: every per-sample FASTQ (R1 and R2 if paired)
    if (alreadyDemuxed) {
        def dir = file(params.fastq_dir)
        if (!dir.isDirectory()) fail("fastq_dir is not a directory: ${params.fastq_dir}")
        def notFound = []
        rows.each { r ->
            def reads = params.paired_end ? ["${r.sample_name}_R1", "${r.sample_name}_R2"] : ["${r.sample_name}"]
            reads.eachWithIndex { base, i ->
                def fq = ['.fastq.gz', '.fastq'].collect { ext -> dir.resolve("${base}${ext}".toString()) }.find { it.exists() }
                if (!fq) {
                    notFound << "${base}.fastq.gz (or .fastq)"
                } else {
                    sampleFastqs << fq
                    if (i == 0) r1Name[r.sample_name] = fq.name
                }
            }
        }
        if (notFound) {
            fail("no FASTQ in ${params.fastq_dir} for these samples in ${params.samplesheet} " +
                 "(file names must be <sample_name>.fastq.gz, or _R1/_R2 when paired_end):\n  " + notFound.join('\n  '))
        }
    } else {
        sampleNames.each { s -> r1Name[s] = params.paired_end ? "${s}_R1.fastq.gz" : "${s}.fastq.gz" }
    }

    // ---- Libraries: which samples are counted against which library ---------
    def libFiles = files(params.library).findAll { it.exists() }
    if (!libFiles) fail("library: ${params.library} matched no files.")
    def libByName = [:]   // library name (file name without .csv) -> file
    libFiles.each { f ->
        def name = f.name.replaceFirst(/\.csv$/, '')
        if (libByName.containsKey(name)) fail("two library files share the name '${name}': ${libByName[name]} and ${f}")
        libByName[name] = f
    }
    def libraryGroups = []   // [name, file, samples] — one mageck count each
    if ('library' in columns) {
        // Each sample is counted only against the library its row names, so
        // every count table (and its QC plots) holds just that library's samples.
        def noLibrary = rows.findAll { !it.library }.collect { it.sample_name }
        if (noLibrary) fail("${params.samplesheet} has a library column, but these samples have no value in it: ${noLibrary.join(', ')}")
        def unknown = rows.findAll { !libByName.containsKey(it.library) }.collect { "${it.sample_name} -> ${it.library}" }
        if (unknown) {
            fail("these samples name a library that isn't among the files matched by library: ${params.library}\n  " +
                 unknown.join('\n  ') + "\nAvailable library names (file name without .csv): ${libByName.keySet().join(', ')}")
        }
        rows.collect { it.library }.unique().each { name ->
            libraryGroups << [name: name, file: libByName[name], samples: rows.findAll { it.library == name }.collect { it.sample_name }]
        }
        def unused = libByName.keySet() - libraryGroups.collect { it.name }
        if (unused) log.warn "library files not assigned to any sample in ${params.samplesheet} (skipped): ${unused.join(', ')}"
    } else {
        // No library column: every sample is counted against every library
        // (the legacy step_2_count.py behaviour; identical for one library).
        libByName.each { name, f -> libraryGroups << [name: name, file: f, samples: sampleNames] }
    }

    // Every sample counted against every library (no `library` column, >1
    // library): a sample with 0 reads in one library is then expected, so
    // MAGECK_COUNT leaves it out of that library's table instead of failing.
    def allVsAll = !('library' in columns) && libraryGroups.size() > 1

    log.info """
    CRISPR_SCREEN  ~  ${params.project_name}
    -----------------------------------------
      samples        : ${sampleNames.size()}
      input          : ${alreadyDemuxed ? "already demultiplexed (fastq_dir: ${hostInput('CRISPR_FASTQ_DIR_HOST', params.fastq_dir)})" : "demultiplex ${hostInput('CRISPR_FASTQ_R1_HOST', params.fastq_r1)}"}
      mode           : ${params.paired_end ? 'paired-end' : 'single-end'}
      libraries      : ${libraryGroups.collect { "${it.name} (${it.samples.size()} samples)" }.join(', ')}
      counting       : ${allVsAll ? 'every sample against every library (no samplesheet library column)' : (libraryGroups.size() > 1 ? "each sample against its samplesheet library" : 'one library')}
      contrasts      : ${params.contrasts ? params.contrasts : '(none — stopping after counting/QC)'}
      outdir         : ${params.outdir}
    """

    return [scripts: stepScripts(), sampleNames: sampleNames, adapterArgs: adapterArgs, barcodeLen: barcodeLen,
            barcodes: alreadyDemuxed ? null : rows.collectEntries { [(it.sample_name): it.barcode] },
            alreadyDemuxed: alreadyDemuxed, sampleFastqs: sampleFastqs, r1Name: r1Name,
            libraryGroups: libraryGroups, allVsAll: allVsAll, images: processContainers()]
}

workflow {

    def cfg = prepareRun()

    // Each step's script: bin/<name>, or this project's replacement (see stepScripts()).
    def sc = cfg.scripts
    def theme = Channel.value(sc.my_theme)

    // Final outputs of the steps nothing else reads (reports, plots): METHODS
    // waits on these and on its own inputs, so it runs after every other step.
    finished = Channel.empty()

    // ---- Step 1: per-sample FASTQs (demultiplex, or take them as given) ----
    if (cfg.alreadyDemuxed) {
        // fastq_dir mode: no demultiplexing, no demux QC.
        sample_fastqs = Channel.fromList(cfg.sampleFastqs)
        demux_log     = Channel.empty()
    } else {
        fastq_r1_ch = Channel.fromPath(params.fastq_r1, checkIfExists: true)
        if (params.paired_end) {
            fastq_r2_ch = Channel.fromPath(params.fastq_r2, checkIfExists: true)
            demux_out = DEMULTIPLEX_PE(fastq_r1_ch, fastq_r2_ch, cfg.adapterArgs, cfg.barcodeLen, sc.demultiplex, sc.barcode_position_qc)
        } else {
            demux_out = DEMULTIPLEX_SE(fastq_r1_ch, cfg.adapterArgs, cfg.barcodeLen, sc.demultiplex, sc.barcode_position_qc)
        }
        sample_fastqs = demux_out.fastqs.flatten()
        demux_log     = demux_out.log

        samplesheet_ch = Channel.fromPath(params.samplesheet)
        demux_stats    = CUTADAPT_STATS(demux_log, samplesheet_ch, sc.cutadapt_stats)
        position_plots = PLOT_BARCODE_POSITION(demux_out.position_csv, sc.plot_barcode_position, theme)
        finished = finished.mix(demux_stats.json, position_plots.plots)
    }

    // ---- QC: FastQC per sample FASTQ, MultiQC summary ----------------------
    fastqc_out = FASTQC(sample_fastqs, sc.fastqc)
    multiqc_out = MULTIQC(fastqc_out.zip.mix(demux_log).collect(), sc.multiqc)
    finished = finished.mix(multiqc_out.report)

    // ---- Step 2: mageck count, one job per library ---------------------------
    library_ch     = Channel.fromList(cfg.libraryGroups.collect { g -> tuple(g.name, g.file) })
    mageck_lib_out = BUILD_MAGECK_LIBRARY(library_ch, sc.build_mageck_library)

    // Each library's own samples, in samplesheet order. FASTQ names come from
    // the samplesheet (not channel emission order), so --fastq and
    // --sample-label stay positionally aligned.
    library_samples_ch = Channel.fromList(cfg.libraryGroups.collect { g ->
        tuple(g.name, g.samples.collect { cfg.r1Name[it] }, g.samples.join(','))
    })

    // ---- Guide position: automatic anchor discovery, or a fixed trim5 -------
    // trim5 unset (the default) -> GUIDE_ANCHOR_DISCOVERY finds this library's
    // constant vector sequence from one of its own samples (no prior
    // knowledge of its sequence or length needed — bin/discover_anchor.py),
    // then GUIDE_ANCHOR_TRIM positions every sample's guide at read position
    // 0 regardless of any variable-length spacing in the amplicon before or
    // after that anchor. Robust to any plasmid design, no per-project tuning.
    // Set trim5 explicitly in config.yaml to bypass both and use the old
    // fixed-offset behavior instead (e.g. to reproduce a specific prior run
    // byte-for-byte).
    def autoAnchor = !params.trim5
    if (autoAnchor) {
        // One representative (first) sample's FASTQ per library, for discovery.
        // + this library's MAGeCK-format guides, so discovery can refine the
        // anchor's 3' end against them (bin/discover_anchor.py).
        discovery_in = sample_fastqs.collect().flatMap { all_fastqs ->
            cfg.libraryGroups.collect { g ->
                def wanted = cfg.r1Name[g.samples[0]]
                tuple(g.name, all_fastqs.find { it.name == wanted })
            }
        }.join(mageck_lib_out.library, failOnMismatch: true, failOnDuplicate: true)
        anchor_out = GUIDE_ANCHOR_DISCOVERY(discovery_in, sc.discover_anchor)

        trim_in  = anchor_out.anchor.join(library_samples_ch, failOnMismatch: true, failOnDuplicate: true)
        trim_out = GUIDE_ANCHOR_TRIM(sample_fastqs.collect(), trim_in, sc.guide_anchor_trim)
        // trim_out.trimmed_fastqs: tuple(lib_name, [<sample>_trimmed.fastq.gz, ...])

        // Each library's mageck count gets only ITS OWN trimmed FASTQs: when
        // a sample is counted against several libraries (no samplesheet
        // `library` column), every library's GUIDE_ANCHOR_TRIM writes the same
        // <sample>_trimmed.fastq.gz name, so pooling them into one input list
        // makes MAGECK_COUNT fail with an input file name collision.
        count_fastqs_ch  = trim_out.trimmed_fastqs
        count_samples_ch = Channel.fromList(cfg.libraryGroups.collect { g ->
            tuple(g.name, g.samples.collect { "${it}_trimmed.fastq.gz" }, g.samples.join(','))
        })
        count_trim5 = Channel.value(0)
    } else {
        count_fastqs_ch   = sample_fastqs.collect().flatMap { all_fastqs ->
            cfg.libraryGroups.collect { g -> tuple(g.name, all_fastqs) }
        }
        count_samples_ch  = library_samples_ch
        count_trim5       = Channel.value(params.trim5)
    }

    // tuple(lib_name, mageck_library, fastq_args, sample_label_str, fastqs)
    count_in = mageck_lib_out.library
        .join(count_samples_ch, failOnMismatch: true, failOnDuplicate: true)
        .join(count_fastqs_ch, failOnMismatch: true, failOnDuplicate: true)

    count_out = MAGECK_COUNT(count_in, cfg.allVsAll, count_trim5, sc.mageck_count)
    // results: tuple(lib_name, count_table, count_normalized, mageck_log, unmapped)
    // sample_reads: tuple(lib_name, <prefix>.sample_reads.tsv)

    // ---- Step 2 QC: how each sample's reads split across libraries -----------
    // Only with more than one library. One task over every library: parallel
    // lists (names, count tables, sample_reads, MAGeCK libraries) built from a
    // single joined channel, sorted by library name, so entry i of each list
    // is the same library and the task hash is stable between runs.
    breakdown_tables = Channel.empty()
    if (cfg.libraryGroups.size() > 1) {
        breakdown_in = count_out.results
            .map { lib_name, count_table, _norm, _log, _unmapped -> tuple(lib_name, count_table) }
            .join(count_out.sample_reads, failOnMismatch: true, failOnDuplicate: true)
            .join(mageck_lib_out.library, failOnMismatch: true, failOnDuplicate: true)
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows -> tuple(rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] }, rows.collect { it[3] }) }
        breakdown = LIBRARY_BREAKDOWN(breakdown_in, cfg.sampleNames.join(','), sc.library_breakdown)
        breakdown_tables = breakdown.table.mix(breakdown.wide)
    }

    // ---- Step 2 QC: guide representation, per library --------------------------
    // Keyed by library name so each count table is paired with its own log and library.
    wrangle_in = count_out.results
        .join(mageck_lib_out.library, failOnMismatch: true, failOnDuplicate: true)
        .map { lib_name, count_table, _norm, mageck_log, unmapped, mageck_library -> tuple(lib_name, count_table, mageck_log, unmapped, mageck_library) }
    guide_tables = GUIDE_REPRESENTATION_WRANGLE(wrangle_in, sc.guide_representation_wrangle)

    // Pair each library's tables with its prefix (taken from <prefix>_counts_long.csv)
    guide_tables_by_pool = guide_tables.tables.map { files ->
        def anchor = files.find { it.name.endsWith('_counts_long.csv') }
        tuple(anchor.name - '_counts_long.csv', files)
    }
    guide_plots = GUIDE_REPRESENTATION_PLOT(guide_tables_by_pool, sc.guide_representation_plot, theme)

    // Pair each library's plots with its prefix (taken from
    // <prefix>_guide_representation.svg — like _ecdf.svg, a plot every run
    // always produces, but unlike _ecdf.svg, no other filename's suffix
    // overlaps it: a library with control guides also produces
    // <prefix>_controls_ecdf.svg, which *also* ends in "_ecdf.svg" — a
    // real bug found this way: files.find() then non-deterministically
    // picked whichever of the two the glob happened to list first,
    // sometimes deriving the wrong prefix from "..._controls_ecdf.svg").
    guide_plots_by_pool = guide_plots.plots.map { files ->
        def anchor = files.find { it.name.endsWith('_guide_representation.svg') }
        tuple(anchor.name - '_guide_representation.svg', files)
    }
    // join(): tuple(prefix, tables) + tuple(prefix, plots) -> tuple(prefix, tables, plots)
    report_in = guide_tables_by_pool.join(guide_plots_by_pool, failOnMismatch: true, failOnDuplicate: true)
    qc_reports = GUIDE_REPRESENTATION_REPORT(report_in, sc.guide_representation_report)
    finished = finished.mix(qc_reports.report)

    // With more than one library, also build one report comparing every
    // pool side by side (facet_wrap per pool + combined headline-metric bar
    // charts) — with a single library this would just repeat its own report.
    if (cfg.libraryGroups.size() > 1) {
        // + the library breakdown tables, which add a "reads per library" section.
        compare_report = GUIDE_REPRESENTATION_COMPARE_REPORT(guide_tables.tables.mix(breakdown_tables).collect(), sc.guide_representation_compare_report, theme)
        finished = finished.mix(compare_report.report)
    }

    // ---- Step 3 (optional): mageck test, one task per contrast -----------------
    // VALIDATE_CONTRASTS sees every library's count table and resolves each
    // contrast to the one table holding both its treatment and control.
    test_evidence = Channel.empty()
    if (params.contrasts) {
        count_tables   = count_out.results.map { lib_name, count_table, _norm, _log, _unmapped -> count_table }
        contrasts_ch   = Channel.fromPath(params.contrasts, checkIfExists: true)
        validated      = VALIDATE_CONTRASTS(count_tables.collect(), contrasts_ch, sc.validate_contrasts)

        tables_by_name = count_tables.map { t -> tuple(t.name, t) }
        // quote: '"' matters — a pattern matching several replicates yields a
        // comma-joined, quoted column list ("A_rep1,A_rep2"); without it the
        // list splits and every later field shifts by one.
        test_in = validated.contrasts
            .splitCsv(header: true, quote: '"')
            .map { row -> tuple(row.count_table, row.name, row.treat_cols, row.control_cols) }
            .combine(tables_by_name, by: 0)
        test_out = MAGECK_TEST(test_in, sc.mageck_test)
        // The resolved contrasts + each contrast's MAGeCK log (<name>.log:
        // MAGeCK's exact command line, version, normalization), for METHODS.md.
        test_evidence = validated.contrasts
            .mix(test_out.results.flatten().filter { it.name.endsWith('.log') })
    } else {
        log.info "params.contrasts not set — stopping after counting + QC (this is normal for a library-representation/pilot run)."
    }

    // ---- Materials and methods: METHODS.md, written last ------------------------
    // bin/write_methods.py describes only the steps this run performed, with
    // the parameters and versions they used, from (1) the run's resolved
    // configuration (runInfo(): inputs, libraries and sample assignment,
    // anchor or fixed trim5, replaced scripts, extra arguments, containers)
    // and (2) what the steps recorded: the tools' own logs (cutadapt and
    // MAGeCK log their exact command line and version, so a replaced script
    // or `<script>_args` is described as it actually ran), the anchors and
    // their discovery/refinement tables, the MAGeCK-format libraries, each
    // library's per-sample read table, the resolved contrasts, and
    // SOFTWARE_VERSIONS' probe of every container image. File names are
    // unique across all of these (prefixed by project/library/contrast), so
    // they are staged together in one evidence/ directory.
    versions = SOFTWARE_VERSIONS(Channel.fromList(cfg.images.values().toSet().sort()), sc.software_versions)

    evidence = demux_log
        .mix(mageck_lib_out.library.map { _lib, lib_file -> lib_file })
        .mix(count_out.results.map { _lib, _table, _norm, mageck_log, _unmapped -> mageck_log })
        .mix(count_out.sample_reads.map { _lib, reads_tsv -> reads_tsv })
        .mix(versions.versions)
        .mix(test_evidence)
    if (autoAnchor) {
        // One sample's anchor-trimming cutadapt log per library (every sample
        // of a library is trimmed with the same command), renamed
        // <library>.anchor_trim.log: sample names repeat across libraries.
        trim_logs = trim_out.logs
            .map { lib_name, logs -> tuple(lib_name, [logs].flatten().sort { it.name }[0]) }
            .collectFile { lib_name, trim_log -> ["${lib_name}.anchor_trim.log", trim_log.text] }
        evidence = evidence
            .mix(anchor_out.anchor.map { _lib, anchor_file -> anchor_file })
            .mix(anchor_out.qc, anchor_out.refinement, trim_logs)
    }

    run_info = Channel.of(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(runInfo(cfg, autoAnchor))))
        .collectFile(name: 'run_info.json', newLine: true)
    METHODS(run_info, evidence.collect(), finished.collect().map { it.size() }, sc.discover_anchor, sc.write_methods)
}
// workflow.onComplete lives in nextflow.config, not here — a nested
// `workflow.onComplete {}` inside this `workflow {}` block resolves
// `workflow` to null at runtime (NullPointerException on .success), and a
// top-level one outside the block is rejected by the strict script parser
// ("Statements cannot be mixed with script declarations").
