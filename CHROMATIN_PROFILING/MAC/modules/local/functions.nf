/*
 * Shared by the modules. Every step runs a script that main.nf resolves from
 * bin/, or from the project's replacement (config.yaml key = script name).
 */

// Interpreter chosen by extension, so a replacement may use another language.
def runner(script) {
    def n = script.name.toLowerCase()
    return n.endsWith('.py') ? 'python3' : n.endsWith('.r') ? 'Rscript' : 'bash'
}

// 'true' / 'false' for a script flag: Boolean true, or true/1/yes in any case, is 'true'; anything else 'false'.
def boolStr(v) {
    def s = v?.toString()?.trim()?.toLowerCase()
    return (v == true || s in ['true', '1', 'yes']) ? 'true' : 'false'
}

// csaw readParam(max.frag): csaw_max_frag, else bowtie2_maxins, else 500. The one rule of every consumer
// (BIN_COUNTS, CSAW, CSAW_CUTOFF_SCAN, CLASSIFY_REGIONS).
def maxFrag() {
    return params.csaw_max_frag ?: params.bowtie2_maxins ?: 500
}

// Histone mark name for comparison: lower case, without a leading H3 / H4 / H2A / H2B / H1 when a K follows
// (H3K27me3 and K27me3 both give k27me3); anything else (H3F3A, IgG) is only lower-cased.
def markKey(antibody) {
    def s = (antibody == null ? '' : antibody.toString().trim().toLowerCase())
    def m = (s =~ /^h(?:3|4|2a|2b|1)(k.*)$/)
    return m.matches() ? m.group(1) : s
}

// True when the samplesheet antibody is in the comma list of broad marks (params.broad_marks), case-insensitively and
// with or without the histone prefix. A blank antibody or an empty list is never broad.
def isBroadMark(antibody, broadMarksCsv) {
    def key = markKey(antibody)
    if (!key || broadMarksCsv == null) return false
    return broadMarksCsv.toString().split(',').collect { m -> markKey(m) }.findAll { m -> m }.contains(key)
}
