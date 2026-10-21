/*
 * Shared by the modules. Every step runs a script that main.nf resolves from
 * bin/, or from the project's replacement (config.yaml key = script name).
 */

// Interpreter chosen by extension, so a replacement may use another language.
def runner(script) {
    def n = script.name.toLowerCase()
    return n.endsWith('.py') ? 'python3' : n.endsWith('.r') ? 'Rscript' : 'bash'
}
