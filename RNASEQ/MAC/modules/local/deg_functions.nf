/*
 * Shared by the deg_*.nf modules (the DEG addition).
 */
include { boolStr } from './functions'

// Single-quoted shell word: safe for spaces, '#', '*', '~', '$' and quotes in values such as a design formula or a plot title.
def shq(v) {
    def s = v == null ? '' : v.toString()
    return "'" + s.replace("'", "'\''") + "'"
}

// 'TRUE' / 'FALSE' for an R script flag (true / 1 / yes in any case is TRUE).
def rbool(v) {
    return boolStr(v).toUpperCase()
}

// `--flag='file'` for an optional file input; empty when the input was given as an empty list (nothing staged).
def optArg(String flag, p) {
    def staged = [p].flatten().findAll { x -> x != null }
    return staged ? "--${flag}=${shq(staged[0].toString())}" : ''
}
