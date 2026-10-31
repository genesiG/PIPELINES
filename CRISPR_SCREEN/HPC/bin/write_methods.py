#!/usr/bin/env python3
"""
write_methods.py
================
Writes METHODS.md: the Materials and Methods for one run of this pipeline,
covering only the steps the run performed, with the parameters and software
versions they actually used.

Nothing is described from what a step does by default. Every value comes from
the run itself:
  - run_info.json (main.nf runInfo()): the resolved configuration — input mode,
    samples, libraries and which samples were counted against each, how the
    guide was located (anchor discovery or a fixed trim5), contrasts, which step
    scripts the project replaced (path + MD5), extra `<script>_args`, the
    container image of every step, and the pipeline and Nextflow versions.
  - the evidence directory (what the steps recorded): cutadapt's and MAGeCK's
    own logs, which hold each tool's exact command line and version (so a flag
    changed by a replaced script or by `<script>_args` is described as it ran),
    the anchors with their discovery and refinement tables, the MAGeCK-format
    libraries, each library's per-sample read table, the resolved contrasts and
    the tool versions SOFTWARE_VERSIONS probed inside each container image.
    Data-dependent behaviour is reported too: samples left out of a library's
    count table, MAGeCK replacing median by total-count normalization, the
    offsets MAGeCK chose with --trim-5 AUTO, duplicated sgRNA sequences, and
    which samples MAGeCK used to estimate variance in each comparison.

Each step's wording is written for this pipeline's default script. A step
whose script the project replaced is still described from its tool's logged
command line where there is one, and is marked "[To complete: ...]" where
only the replacement's author knows what it does.

Keep this file in step with bin/: when a step's behaviour, defaults or outputs
change, update its section here.

A section that fails to generate is replaced by a visible note (and the error
is printed to stderr) rather than failing the run at its very last step.

USAGE:
  python3 write_methods.py --run-info=run_info.json --evidence-dir=evidence \
      [--discover-script=discover_anchor.py] [--output=METHODS.md] \
      [--run-info-out=run_info.json] [--versions-out=software_versions.tsv]

  --discover-script  the discover_anchor script this run used; its argparse
                     defaults are read so the quoted parameters are the ones
                     that applied (plus any discover_anchor_args)

Evidence file names (as main.nf stages them):
  cutadapt.log                              demultiplexing (cutadapt) log
  <library>_mageck_library.txt              MAGeCK-format library
  <prefix>_mageck_count.log                 mageck count log (<prefix> = <project>_<library>)
  <prefix>.sample_reads.tsv                 sample, fastq, reads, mapped, status
  <project>_<library>_anchor.txt            anchor used for the library
  <project>_<library>_anchor_discovery.csv  discovery trace
  <project>_<library>_anchor_refinement.csv refinement candidates
  <library>.anchor_trim.log                 one sample's anchor-trimming cutadapt log
  validated_contrasts.csv                   resolved contrasts
  <contrast name>.log                       mageck test log
  *.versions.tsv                            image, tool, version

OUTPUTS:
  METHODS.md              Materials and Methods, software table, references,
                          exact commands and input checksums
  run_info.json           copy of the run's resolved configuration
  software_versions.tsv   image, tool, version (all images probed)
"""

import argparse
import csv
import json
import math
import os
import re
import sys
import traceback
from datetime import datetime

# Same rule as guide_representation_wrangle.R's is_control.
CONTROL_PATTERN = re.compile(r"Control|NTC|Non-targeting", re.IGNORECASE)

# ---------------------------------------------------------------------------
# References (author-year). Only those cited end up in METHODS.md.
# ---------------------------------------------------------------------------
REFERENCES = {
    "andrews2010": ("Andrews, 2010",
                    "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. "
                    "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"),
    "benjamini1995": ("Benjamini & Hochberg, 1995",
                      "Benjamini Y, Hochberg Y. (1995). Controlling the false discovery rate: a practical and "
                      "powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol 57:289–300."),
    "daveiga2017": ("da Veiga Leprevost et al., 2017",
                    "da Veiga Leprevost F, Grüning BA, Alves Aflitos S, et al. (2017). BioContainers: an "
                    "open-source and community-driven framework for software standardization. "
                    "Bioinformatics 33:2580–2582. doi:10.1093/bioinformatics/btx192"),
    "ditommaso2017": ("Di Tommaso et al., 2017",
                      "Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. (2017). Nextflow "
                      "enables reproducible computational workflows. Nat Biotechnol 35:316–319. "
                      "doi:10.1038/nbt.3820"),
    "ewels2016": ("Ewels et al., 2016",
                  "Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for "
                  "multiple tools and samples in a single report. Bioinformatics 32:3047–3048. "
                  "doi:10.1093/bioinformatics/btw354"),
    "kurtzer2017": ("Kurtzer et al., 2017",
                    "Kurtzer GM, Sochat V, Bauer MW. (2017). Singularity: scientific containers for mobility of "
                    "compute. PLoS One 12:e0177459. doi:10.1371/journal.pone.0177459"),
    "li2014": ("Li et al., 2014",
               "Li W, Xu H, Xiao T, Cong L, Love MI, Zhang F, Irizarry RA, Liu JS, Brown M, Liu XS. (2014). "
               "MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 "
               "knockout screens. Genome Biol 15:554. doi:10.1186/s13059-014-0554-4"),
    "li2015": ("Li et al., 2015",
               "Li W, Köster J, Xu H, Chen CH, Xiao T, Liu JS, Brown M, Liu XS. (2015). Quality control, "
               "modeling, and visualization of CRISPR screens with MAGeCK-VISPR. Genome Biol 16:281. "
               "doi:10.1186/s13059-015-0843-6"),
    "martin2011": ("Martin, 2011",
                   "Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. "
                   "EMBnet.journal 17:10–12. doi:10.14806/ej.17.1.200"),
    "rcore": ("R Core Team, {year}",
              "R Core Team ({year}). R: a language and environment for statistical computing. R Foundation for "
              "Statistical Computing, Vienna, Austria. https://www.R-project.org/"),
    "wickham2016": ("Wickham, 2016",
                    "Wickham H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York. "
                    "ISBN 978-3-319-24277-4. https://ggplot2.tidyverse.org"),
    "wickham2019": ("Wickham et al., 2019",
                    "Wickham H, Averick M, Bryan J, et al. (2019). Welcome to the tidyverse. J Open Source Softw "
                    "4:1686. doi:10.21105/joss.01686"),
}

# Which scripts each process runs (main.nf / modules), and a label for each.
PROCESS_SCRIPTS = {
    "DEMULTIPLEX_SE": ["demultiplex", "barcode_position_qc"],
    "DEMULTIPLEX_PE": ["demultiplex", "barcode_position_qc"],
    "CUTADAPT_STATS": ["cutadapt_stats"],
    "PLOT_BARCODE_POSITION": ["plot_barcode_position", "my_theme"],
    "FASTQC": ["fastqc"],
    "MULTIQC": ["multiqc"],
    "BUILD_MAGECK_LIBRARY": ["build_mageck_library"],
    "GUIDE_ANCHOR_DISCOVERY": ["discover_anchor"],
    "GUIDE_ANCHOR_TRIM": ["guide_anchor_trim"],
    "MAGECK_COUNT": ["mageck_count"],
    "LIBRARY_BREAKDOWN": ["library_breakdown"],
    "GUIDE_REPRESENTATION_WRANGLE": ["guide_representation_wrangle"],
    "GUIDE_REPRESENTATION_PLOT": ["guide_representation_plot", "my_theme"],
    "GUIDE_REPRESENTATION_REPORT": ["guide_representation_report"],
    "GUIDE_REPRESENTATION_COMPARE_REPORT": ["guide_representation_compare_report", "my_theme"],
    "VALIDATE_CONTRASTS": ["validate_contrasts"],
    "MAGECK_TEST": ["mageck_test"],
    "SOFTWARE_VERSIONS": ["software_versions"],
    "METHODS": ["write_methods"],
}
SCRIPT_LABELS = {
    "demultiplex": "Sample demultiplexing",
    "barcode_position_qc": "Barcode-position QC table",
    "cutadapt_stats": "Demultiplexing statistics",
    "plot_barcode_position": "Barcode-position plot",
    "my_theme": "Plot theme",
    "fastqc": "Read quality (FastQC)",
    "multiqc": "QC summary (MultiQC)",
    "build_mageck_library": "Library conversion",
    "discover_anchor": "Anchor discovery",
    "guide_anchor_trim": "Anchor trimming",
    "mageck_count": "Read counting",
    "library_breakdown": "Reads per library",
    "guide_representation_wrangle": "Library representation QC (metrics)",
    "guide_representation_plot": "Library representation QC (plots)",
    "guide_representation_report": "Library representation QC (report)",
    "guide_representation_compare_report": "Library comparison report",
    "validate_contrasts": "Contrast resolution",
    "mageck_test": "Statistical testing",
    "software_versions": "Software versions",
    "write_methods": "This METHODS.md",
}
STEP_TOOLS = {   # tool -> [(process, what it was used for)]
    "cutadapt": [("DEMULTIPLEX_SE", "sample demultiplexing"), ("DEMULTIPLEX_PE", "sample demultiplexing"),
                 ("GUIDE_ANCHOR_TRIM", "anchor trimming")],
    "MAGeCK": [("MAGECK_COUNT", "read counting"), ("MAGECK_TEST", "gene-level testing")],
    "FastQC": [("FASTQC", "read quality assessment")],
    "MultiQC": [("MULTIQC", "QC report summary")],
}
PYTHON_PROCESSES = ["BUILD_MAGECK_LIBRARY", "GUIDE_ANCHOR_DISCOVERY", "LIBRARY_BREAKDOWN", "VALIDATE_CONTRASTS",
                    "CUTADAPT_STATS", "METHODS"]
R_PROCESSES = ["GUIDE_REPRESENTATION_WRANGLE", "GUIDE_REPRESENTATION_PLOT", "GUIDE_REPRESENTATION_REPORT",
               "GUIDE_REPRESENTATION_COMPARE_REPORT", "PLOT_BARCODE_POSITION"]
R_PACKAGES = ["dplyr", "tidyr", "ggplot2", "scales"]   # = bin/software_versions.sh R_PACKAGES

# discover_anchor.py's extension starts at k = 4 (range(4, max_length + 1)).
DISCOVERY_START_K = 4


# ---------------------------------------------------------------------------
# Small text helpers
# ---------------------------------------------------------------------------
def oxford(items) -> str:
    items = [str(x) for x in items]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def count_of(n: int, singular: str, plural: str = "") -> str:
    return f"{n:,} {singular if n == 1 else (plural or singular + 's')}"


def code(text) -> str:
    text = str(text)
    fence = "``" if "`" in text else "`"
    return f"{fence}{text}{fence}"


def num(value) -> str:
    """0.1 -> '0.1', 8.0 -> '8', 200000 -> '200,000'."""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return f"{value:,}" if isinstance(value, int) else str(value)


def natural_key(name: str):
    """Pool_2 before Pool_10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(name))]


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def table(headers, rows) -> str:
    def cell(c):
        return str(c).replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(cell(h) for h in headers) + " |",
             "|" + "|".join("---" for _ in headers) + "|"]
    lines += ["| " + " | ".join(cell(c) for c in r) + " |" for r in rows]
    return "\n".join(lines)


def to_complete(what: str) -> str:
    return f"**[To complete: {what}]**"


def read_tsv(path, delimiter="\t"):
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


# ---------------------------------------------------------------------------
# Command-line parsing (the tools log their argv as space-joined tokens)
# ---------------------------------------------------------------------------
def parse_argv(argv, value_opts, multi_opts=()):
    """-> (options [(name, value|list|None)], positional). Handles --opt=value
    and attached short values (-e0)."""
    opts, positional = [], []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok.startswith("--") and "=" in tok:
            name, value = tok.split("=", 1)
            opts.append((name, value))
            i += 1
        elif tok.startswith("-") and len(tok) > 1 and to_float(tok) is None:
            if tok in multi_opts:
                values, i = [], i + 1
                while i < len(argv) and not (argv[i].startswith("-") and to_float(argv[i]) is None):
                    values.append(argv[i])
                    i += 1
                opts.append((tok, values))
            elif not tok.startswith("--") and len(tok) > 2 and tok[:2] in value_opts:
                opts.append((tok[:2], tok[2:]))
                i += 1
            elif tok in value_opts and i + 1 < len(argv):
                opts.append((tok, argv[i + 1]))
                i += 2
            else:
                opts.append((tok, None))
                i += 1
        else:
            positional.append(tok)
            i += 1
    return opts, positional


def last(opts, *names):
    found = [v for n, v in opts if n in names]
    return found[-1] if found else None


def present(opts, *names) -> bool:
    return any(n in names for n, _ in opts)


def show(opts, names) -> str:
    """The options among `names` as they were given, e.g. '-e 0 --no-indels'."""
    parts = []
    for n, v in opts:
        if n in names:
            parts.append(n if v is None else f"{n} {' '.join(v) if isinstance(v, list) else v}")
    return " ".join(parts)


# ---- cutadapt ---------------------------------------------------------------
CUTADAPT_ADAPTER_OPTS = {"-a", "-g", "-b", "-A", "-G", "-B", "--adapter", "--front", "--anywhere"}
CUTADAPT_MATCH_OPTS = {"-e", "--error-rate", "--errors", "-O", "--overlap", "--no-indels",
                       "--discard-untrimmed", "--trimmed-only"}
CUTADAPT_HOUSEKEEPING = {"-o", "--output", "-p", "--paired-output", "--info-file", "-j", "--cores", "--json",
                         "--report", "--quiet", "--debug", "-Z", "--compression-level", "--buffer-size"}
CUTADAPT_QUALITY = {   # option -> description ({v} = its value)
    "-q": "3′ quality trimming at Phred {v}", "--quality-cutoff": "3′ quality trimming at Phred {v}",
    "--nextseq-trim": "two-colour (NextSeq) quality trimming at Phred {v}",
    "-m": "removal of reads shorter than {v} nt", "--minimum-length": "removal of reads shorter than {v} nt",
    "-M": "removal of reads longer than {v} nt", "--maximum-length": "removal of reads longer than {v} nt",
    "--max-n": "removal of reads with more than {v} N bases",
    "--max-ee": "removal of reads with more than {v} expected errors",
    "--max-expected-errors": "removal of reads with more than {v} expected errors",
    "--max-average-error-rate": "removal of reads with an average error rate above {v}",
    "--max-aer": "removal of reads with an average error rate above {v}",
    "--trim-n": "trimming of N bases at read ends",
    "-u": "removal of a fixed {v} bases", "--cut": "removal of a fixed {v} bases",
    "-l": "shortening of reads to {v} nt", "--length": "shortening of reads to {v} nt",
    "--poly-a": "poly-A tail trimming",
    "--discard-casava": "removal of reads that failed the Illumina chastity filter",
}
CUTADAPT_VALUE_OPTS = (CUTADAPT_ADAPTER_OPTS | {"-e", "--error-rate", "--errors", "-O", "--overlap", "-n", "--times",
                       "--action", "-o", "--output", "-p", "--paired-output", "--info-file", "-j", "--cores",
                       "--json", "--report", "-q", "--quality-cutoff", "--nextseq-trim", "--quality-base", "-m",
                       "--minimum-length", "-M", "--maximum-length", "--max-n", "--max-ee",
                       "--max-expected-errors", "--max-average-error-rate", "--max-aer", "-u", "--cut", "-U",
                       "-l", "--length", "--length-tag", "--strip-suffix", "-x", "--prefix", "-y", "--suffix",
                       "--rename", "-r", "--rest-file", "--wildcard-file", "--too-short-output",
                       "--too-long-output", "--untrimmed-output", "--untrimmed-paired-output",
                       "--too-short-paired-output", "--too-long-paired-output", "--pair-filter",
                       "--compression-level", "--buffer-size"})


class CutadaptRun:
    """One cutadapt invocation, from its log ('This is cutadapt X' +
    'Command line parameters: ...')."""

    def __init__(self, text: str):
        self.version = self.argv = None
        m = re.search(r"^This is cutadapt (\S+)", text or "", re.M)
        self.version = m.group(1) if m else None
        m = re.search(r"^Command line parameters: (.*)$", text or "", re.M)
        if not m:
            return
        self.argv = m.group(1).split()
        self.opts, self.positional = parse_argv(self.argv, CUTADAPT_VALUE_OPTS)
        self.adapters = [(n, v) for n, v in self.opts if n in CUTADAPT_ADAPTER_OPTS]
        e = last(self.opts, "-e", "--error-rate", "--errors")
        self.error_rate = to_float(e)
        o = last(self.opts, "-O", "--overlap")
        self.overlap = int(o) if o is not None and o.isdigit() else None
        self.no_indels = present(self.opts, "--no-indels")
        self.discard_untrimmed = present(self.opts, "--discard-untrimmed", "--trimmed-only")
        self.paired = present(self.opts, "-p", "--paired-output")
        self.quality = [(n, v) for n, v in self.opts if n in CUTADAPT_QUALITY]
        known = CUTADAPT_ADAPTER_OPTS | CUTADAPT_MATCH_OPTS | CUTADAPT_HOUSEKEEPING | set(CUTADAPT_QUALITY)
        self.other = [(n, v) for n, v in self.opts if n not in known]

    @property
    def ok(self) -> bool:
        return self.argv is not None

    def adapter_seqs(self):
        """[(name, sequence, anchored)] of the 5'/3' adapters, 'name=' and '^'/'$' removed."""
        out = []
        for _, spec in self.adapters:
            name, seq = (spec.split("=", 1) if "=" in spec else (None, spec))
            if seq.startswith("file:"):
                continue
            anchored = seq.startswith("^")
            out.append((name, seq.strip("^$").split(";")[0].upper(), anchored))
        return out

    def quality_phrases(self):
        return [CUTADAPT_QUALITY[n].format(v=v) for n, v in self.quality]

    def other_text(self):
        return " ".join(n if v is None else f"{n} {v}" for n, v in self.other)


def error_phrase(e, no_indels, lengths, flags) -> str:
    """How many errors a match could have, for sequences of the given lengths."""
    kinds = "mismatches" if no_indels else "mismatches or insertions/deletions"
    shown = f" ({code(flags)})" if flags else ""
    if e is None:
        e, shown = 0.1, " (cutadapt's default error rate, 0.1)"
    if e == 0:
        return ("with no mismatches, insertions or deletions" if no_indels else "with no errors") + shown
    if e >= 1:
        return f"with up to {int(e)} {kinds}{shown}"
    per_len = {L: int(math.floor(e * L)) for L in sorted(set(lengths))}
    if len(set(per_len.values())) == 1:
        k = next(iter(per_len.values()))
        span = oxford([f"{L}" for L in per_len]) + "-nt"
        return f"with an error rate of at most {num(e)}, i.e. up to {k} {kinds} in the {span} sequence{shown}"
    detail = "; ".join(f"{k} for {L} nt" for L, k in per_len.items())
    return f"with an error rate of at most {num(e)} (⌊{num(e)} × length⌋ {kinds}: {detail}){shown}"


# ---- MAGeCK -------------------------------------------------------------------
MAGECK_VALUE_OPTS = {"-l", "--list-seq", "-k", "--count-table", "--norm-method", "--control-sgrna", "--control-gene",
                     "--sample-label", "-n", "--output-prefix", "--count-pair", "--trim-5", "--sgrna-len",
                     "--day0-label", "--gmt-file", "-t", "--treatment-id", "-c", "--control-id",
                     "--gene-test-fdr-threshold", "--adjust-method", "--variance-estimation-samples",
                     "--sort-criteria", "--remove-zero", "--remove-zero-threshold", "--gene-lfc-method",
                     "--skip-gene", "--additional-rra-parameters", "--cnv-norm", "--cell-line", "--cnv-est"}
MAGECK_MULTI_OPTS = {"--fastq", "--fastq-2"}
MAGECK_COUNT_DESCRIBED = {"-l", "--list-seq", "--fastq", "--sample-label", "-n", "--output-prefix", "--trim-5",
                          "--sgrna-len", "--norm-method", "--unmapped-to-file"}
MAGECK_TEST_DESCRIBED = {"-k", "--count-table", "-t", "--treatment-id", "-c", "--control-id", "-n", "--output-prefix",
                         "--norm-method", "--adjust-method", "--paired", "--control-sgrna", "--control-gene",
                         "--remove-zero", "--remove-zero-threshold", "--gene-lfc-method", "--gene-test-fdr-threshold",
                         "--variance-estimation-samples", "--sort-criteria"}


class MageckRun:
    """One MAGeCK invocation, from its log ('Parameters: ...', 'Welcome to MAGeCK vX')."""

    def __init__(self, text: str):
        text = text or ""
        self.version = self.subcommand = None
        self.opts, self.positional, self.argv = [], [], None
        m = re.search(r"Welcome to MAGeCK v(\S+?)\.?\s+Command", text)
        self.version = m.group(1) if m else None
        m = re.search(r"Parameters: (.*?)\s*$", text, re.M)
        if m:
            argv = m.group(1).split()
            idx = next((i for i, t in enumerate(argv) if os.path.basename(t) == "mageck"), None)
            if idx is not None and idx + 1 < len(argv):
                self.subcommand = argv[idx + 1]
                self.argv = ["mageck"] + argv[idx + 1:]
                self.opts, self.positional = parse_argv(argv[idx + 2:], MAGECK_VALUE_OPTS, MAGECK_MULTI_OPTS)
        # Per FASTQ file: reads, mapped, AUTO trim result.
        self.files, current = [], None
        for line in text.splitlines():
            m = re.search(r"Parsing FASTQ file (.*?)\.\.\.", line)
            if m:
                current = {"fastq": m.group(1).strip(), "total": None, "mapped": None}
                self.files.append(current)
            m = re.search(r"Determining the trim-5 length of FASTQ file (.*?)\.\.\.", line)
            if m:
                self.files.append({"fastq": m.group(1).strip(), "auto": None, "probe": True})
            m = re.search(r"Auto determination of trim5 results: (\S+)", line)
            if m and self.files:
                self.files[-1]["auto"] = m.group(1)
            m = re.search(r"INFO .*: Total: (\d+)", line)
            if m and current:
                current["total"] = int(m.group(1))
            m = re.search(r"INFO .*: Mapped: (\d+)", line)
            if m and current:
                current["mapped"] = int(m.group(1))
        self.auto_trim = {f["fastq"]: f["auto"] for f in self.files if f.get("probe") and f.get("auto")}
        m = re.search(r"There are (\d+) sgRNAs with duplicated sequences", text)
        self.duplicates = int(m.group(1)) if m else 0
        # MAGeCK falls back from median to total-count normalization for the whole table.
        self.total_norm = []   # [(sample index, reason)]
        for m in re.finditer(r"Sample (\d+) has zero median count, so median normalization is not possible", text):
            self.total_norm.append((int(m.group(1)), "a median count of zero"))
        for m in re.finditer(r"Sample (\d+) has too many zero-count sgRNAs \(([\d.]+)\), and median normalization "
                             r"is unstable\. Switch", text):
            self.total_norm.append((int(m.group(1)), f"{float(m.group(2)) * 100:.0f}% of sgRNAs with zero reads"))

    @property
    def ok(self) -> bool:
        return self.argv is not None

    def opt(self, *names):
        return last(self.opts, *names)

    def labels(self):
        lab = self.opt("--sample-label")
        return lab.split(",") if lab else []

    def other(self, described):
        return [(n, v) for n, v in self.opts if n not in described]

    def short_command(self, max_items=3) -> str:
        """The command with long FASTQ / label lists collapsed."""
        parts = ["mageck", self.subcommand]
        for n, v in self.opts:
            if isinstance(v, list):
                parts.append(f"{n} {' '.join(v)}" if len(v) <= max_items else f"{n} <{len(v)} FASTQ files>")
            elif n == "--sample-label" and v and v.count(",") >= max_items:
                parts.append(f"{n} <{v.count(',') + 1} sample names>")
            else:
                parts.append(n if v is None else f"{n} {v}")
        return " ".join(parts + self.positional)


# ---------------------------------------------------------------------------
# The document
# ---------------------------------------------------------------------------
class Methods:
    def __init__(self, info: dict, evidence_dir: str, discovery_defaults: dict):
        self.info = info
        self.dir = evidence_dir
        self.names = set(os.listdir(evidence_dir)) if os.path.isdir(evidence_dir) else set()
        self.cited, self.ref_fmt = [], {}
        self.n_tables = 0
        self.params = info.get("params") or {}
        self.project = self.params.get("project_name") or "project"
        self.libraries = sorted(info.get("libraries") or [], key=lambda x: natural_key(x["name"]))
        self.lib_names = [lib["name"] for lib in self.libraries]
        self.samples = info.get("samples") or []
        self.input = info.get("input") or {}
        self.demux = self.input.get("mode") == "demultiplex"
        self.paired = bool(self.input.get("paired_end"))
        self.position = (info.get("guide_position") or {}).get("mode", "fixed")
        self.scripts = info.get("scripts") or {}
        self.extra_args = info.get("extra_args") or {}
        self.containers = info.get("containers") or {}
        self.stub = bool((info.get("nextflow") or {}).get("stub_run"))
        self.has_contrasts = bool(info.get("contrasts"))
        self.processes = self._processes()
        self.versions = self._versions()
        self.discover_defaults = discovery_defaults

        # Per-library evidence.
        self.lib_seqs, self.count_runs, self.sample_reads = {}, {}, {}
        self.anchors, self.discovery, self.refinement, self.trim_runs = {}, {}, {}, {}
        for lib in self.lib_names:
            self.lib_seqs[lib] = self._library(f"{lib}_mageck_library.txt")
            prefix = f"{self.project}_{lib}"
            self.count_runs[lib] = MageckRun(self._text(f"{prefix}_mageck_count.log"))
            self.sample_reads[lib] = self._rows(f"{prefix}.sample_reads.tsv")
            anchor = self._text(f"{self.project}_{lib}_anchor.txt")
            self.anchors[lib] = anchor.strip().upper() if anchor and anchor.strip() else None
            self.discovery[lib] = self._rows(f"{self.project}_{lib}_anchor_discovery.csv", ",")
            self.refinement[lib] = self._rows(f"{self.project}_{lib}_anchor_refinement.csv", ",")
            trim = self._text(f"{lib}.anchor_trim.log")
            self.trim_runs[lib] = CutadaptRun(trim) if trim else None
        demux = self._text("cutadapt.log")
        self.demux_run = CutadaptRun(demux) if demux else None
        self.contrasts = self._rows("validated_contrasts.csv", ",")
        self.test_runs = {c["name"]: MageckRun(self._text(f"{c['name']}.log")) for c in self.contrasts}

    # ---- evidence access ------------------------------------------------------
    def _path(self, name):
        return os.path.join(self.dir, name) if name in self.names else None

    def _text(self, name):
        p = self._path(name)
        if not p:
            return None
        with open(p, encoding="utf-8", errors="replace") as fh:
            return fh.read()

    def _rows(self, name, delimiter="\t"):
        p = self._path(name)
        if not p or os.path.getsize(p) == 0:
            return []
        try:
            return read_tsv(p, delimiter)
        except (csv.Error, UnicodeDecodeError):
            return []

    def _library(self, name):
        """[(sgRNA_id, sequence, gene)] of a MAGeCK-format library."""
        rows = []
        p = self._path(name)
        if not p:
            return None
        with open(p, encoding="utf-8-sig", errors="replace") as fh:
            next(fh, None)
            for line in fh:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 3 and parts[1].strip():
                    rows.append((parts[0], parts[1].strip().upper(), parts[2]))
        return rows

    def _versions(self):
        out = {}
        for name in sorted(self.names):
            if name.endswith(".versions.tsv"):
                for r in self._rows(name):
                    if r.get("tool") and r.get("version"):
                        out.setdefault(r.get("image", ""), {})[r["tool"]] = r["version"]
        return out

    def _processes(self):
        p = []
        if self.demux:
            p += ["DEMULTIPLEX_PE" if self.paired else "DEMULTIPLEX_SE", "CUTADAPT_STATS", "PLOT_BARCODE_POSITION"]
        p += ["FASTQC", "MULTIQC", "BUILD_MAGECK_LIBRARY"]
        if self.position == "anchor":
            p += ["GUIDE_ANCHOR_DISCOVERY", "GUIDE_ANCHOR_TRIM"]
        p += ["MAGECK_COUNT", "GUIDE_REPRESENTATION_WRANGLE", "GUIDE_REPRESENTATION_PLOT", "GUIDE_REPRESENTATION_REPORT"]
        if len(self.libraries) > 1:
            p += ["LIBRARY_BREAKDOWN", "GUIDE_REPRESENTATION_COMPARE_REPORT"]
        if self.has_contrasts:
            p += ["VALIDATE_CONTRASTS", "MAGECK_TEST"]
        return p + ["SOFTWARE_VERSIONS", "METHODS"]

    # ---- shared lookups ---------------------------------------------------------
    def cite(self, key, **fmt) -> str:
        """In-text citation; fmt fills placeholders such as R's {year}."""
        if key not in self.cited:
            self.cited.append(key)
        self.ref_fmt.update(fmt)
        return REFERENCES[key][0].format(**self.ref_fmt)

    def table_no(self) -> str:
        self.n_tables += 1
        return f"Table {self.n_tables}"

    def version(self, tool, process):
        image = self.containers.get(process)
        v = self.versions.get(image, {}).get(tool)
        if v:
            return v
        for tools in self.versions.values():   # same tool in another image
            if tool in tools:
                return tools[tool]
        return None

    def v(self, tool, process, logged=None) -> str:
        """' v1.2' — the version the tool logged itself if it did, else the probed one."""
        ver = logged or self.version(tool, process)
        return f" v{ver}" if ver else " (version not recorded)"

    def replaced(self, script) -> bool:
        return bool((self.scripts.get(script) or {}).get("replaced"))

    def replaced_note(self, script) -> str:
        if not self.replaced(script):
            return ""
        return (f" This step ran the project-specific script {code(self.scripts[script]['path'])} instead of the "
                f"pipeline's default (see Software and reproducibility).")

    def extra_note(self, script, tool_logged=True) -> str:
        """Extra arguments of a helper script (a tool's own extras are read from its logged command)."""
        args = self.extra_args.get(script)
        if not args or tool_logged:
            return ""
        return f" Additional arguments: {code(args)}."

    def lib_lengths(self, lib):
        seqs = self.lib_seqs.get(lib) or []
        return sorted({len(s) for _, s, _ in seqs})

    def active_samples(self, lib):
        """Samples in the library's count table with at least one counted read."""
        rows = self.sample_reads.get(lib) or []
        if rows:
            return [r["sample"] for r in rows if r.get("status") == "counted" and (to_float(r.get("mapped")) or 0) > 0]
        return list(next((x["samples"] for x in self.libraries if x["name"] == lib), []))

    def dropped_samples(self, lib, status="dropped_zero_mapped"):
        return [r["sample"] for r in (self.sample_reads.get(lib) or []) if r.get("status") == status]

    def other_drop_statuses(self):
        """Statuses other than counted / dropped_zero_mapped (e.g. from a replaced counting script)."""
        return sorted({r.get("status") for rows in self.sample_reads.values() for r in rows}
                      - {"counted", "dropped_zero_mapped", None, ""})

    # ---- sections ---------------------------------------------------------------
    def header(self):
        nf = self.info.get("nextflow") or {}
        pipe = self.info.get("pipeline") or {}
        start = nf.get("start") or ""
        try:
            when = datetime.fromisoformat(start[:19]).strftime("%Y-%m-%d %H:%M")
        except ValueError:
            when = start
        commit = pipe.get("commit")
        version = f"v{pipe.get('version')}" if pipe.get("version") else ""
        if commit:
            version += f", commit {commit}" + (" with uncommitted changes" if pipe.get("modified") else "")
        config = nf.get("config_file") or "config.yaml"
        blocks = [f"# Materials and Methods — {self.project}"]
        quote = (f"> Generated automatically by the {pipe.get('name') or 'CRISPR_SCREEN'} pipeline ({version}) "
                 f"at the end of run `{nf.get('run_name', '?')}` (started {when}; configuration {code(config)}, "
                 f"results in {code(self.params.get('outdir', './results'))}).\n>\n"
                 f"> It describes the computational analysis only, and only the steps this run performed, with "
                 f"the parameters and software versions they used. Library cloning, cell culture, genomic DNA "
                 f"extraction, sgRNA amplification and sequencing need to be described separately. Text marked "
                 f"**[To complete: …]** needs information the pipeline does not have. Values in `code` are the "
                 f"exact options given to each tool; Appendix A lists the full commands.")
        if self.stub:
            quote += ("\n>\n> **This file comes from a `-stub-run`: every step wrote placeholder outputs, so it does "
                      "not describe a real analysis.**")
        blocks.append(quote)
        blocks.append("## CRISPR screen sequencing analysis")
        return blocks

    def demultiplexing(self):
        n = len(self.samples)
        if not self.demux:
            reads = ("paired-end; read 1 of each pair was used for sgRNA counting" if self.paired
                     else "single-end")
            return ["### Sequencing reads",
                    f"Sequencing reads were provided as already demultiplexed per-sample FASTQ files for "
                    f"{count_of(n, 'sample')} ({reads}), so no demultiplexing was performed in this analysis. "
                    + to_complete("how the reads were demultiplexed into samples (for example by the sequencing "
                                  "facility from index reads, or by inline barcode, and with which tool and "
                                  "mismatch tolerance)")]
        run = self.demux_run
        title = "### Sample demultiplexing"
        if not run or not run.ok:
            return [title, f"Reads were demultiplexed into {count_of(n, 'sample')} by inline barcode."
                    + self.replaced_note("demultiplex") + " "
                    + to_complete("the demultiplexing method; no cutadapt log was found for this run")]
        cut = self.cite("martin2011")
        seqs = run.adapter_seqs()
        lengths = sorted({len(s) for _, s, _ in seqs}) or [0]
        anchored = bool(seqs) and all(a for _, _, a in seqs)
        len_text = f"{oxford(lengths)}-nt" if lengths != [0] else ""
        text = (f"Reads{' (read 1 of each pair)' if run.paired else ''} were assigned to samples "
                f"({count_of(n, 'sample')}) by their inline {len_text} sample barcode with cutadapt"
                f"{self.v('cutadapt', self.processes[0], run.version)} ({cut}) in demultiplexing mode. ")
        if anchored:
            text += ("Each sample's barcode was given as a 5′ adapter anchored at the read start "
                     f"({code('-g <sample>=^<barcode>')}), so a barcode had to occupy the first bases of the read")
        else:
            text += (f"Each sample's barcode was given as a 5′ adapter ({code('-g <sample>=<barcode>')}), which "
                     f"cutadapt searches for anywhere in the read")
        flags = show(run.opts, {"-e", "--error-rate", "--errors", "--no-indels"})
        overlap = self._overlap_phrase(run, lengths, "barcode", anchored)
        text += (f"; a read was assigned to a sample only if it contained {overlap} "
                 f"{error_phrase(run.error_rate, run.no_indels, lengths, flags)}. The barcode and every base 5′ of "
                 f"it were removed from the read. ")
        if run.discard_untrimmed:
            text += f"Reads in which no barcode was found were discarded ({code('--discard-untrimmed')})."
        else:
            text += "Reads in which no barcode was found were set aside and not analysed further."
        if run.paired:
            text += (" Barcodes were searched in read 1 only, and read 2 was assigned to the same sample as its "
                     "mate; only read 1 was used for sgRNA counting.")
        if run.quality:
            text += f" In addition, cutadapt applied {oxford(run.quality_phrases())} ({code(show(run.opts, set(CUTADAPT_QUALITY)))})."
        if run.other:
            text += f" Other cutadapt options: {code(run.other_text())}."
        if len(seqs) != n and seqs:
            text += f" (The cutadapt command listed {count_of(len(seqs), 'barcode')} for {count_of(n, 'sample')}.)"
        return [title, text + self.replaced_note("demultiplex")]

    def _overlap_phrase(self, run, lengths, what, anchored) -> str:
        """What part of the sample's barcode a read had to contain."""
        if anchored:
            return f"that sample's complete {what}"
        o = run.overlap
        if o is None:
            return f"at least 3 nt of that sample's {what} (cutadapt's default minimum overlap)"
        short = sorted({L - o for L in lengths if L > o})
        if not short:
            return f"that sample's complete {what} ({code('--overlap ' + str(o))})"
        return (f"at least {o} nt of that sample's {what} ({code('--overlap ' + str(o))}; a {what} truncated by "
                f"up to {oxford(short)} nt at the start of the read was also accepted)")

    def read_quality(self):
        fq = self.v("FastQC", "FASTQC")
        mq = self.v("MultiQC", "MULTIQC")
        text = (f"Read quality of every per-sample FASTQ file was assessed with FastQC{fq} "
                f"({self.cite('andrews2010')}), and the FastQC reports"
                f"{', together with the cutadapt demultiplexing report,' if self.demux else ''} were summarised "
                f"with MultiQC{mq} ({self.cite('ewels2016')}). These reports were used for inspection only. ")
        steps = []
        if self.demux_run and self.demux_run.ok and self.demux_run.quality:
            steps.append("demultiplexing")
        if any(r and r.ok and r.quality for r in self.trim_runs.values()):
            steps.append("anchor trimming")
        if steps:
            text += (f"Quality-based read processing was limited to the cutadapt options described for "
                     f"{oxford(steps)}.")
        else:
            text += ("No base-quality trimming or quality-based read filtering was applied: reads were kept or "
                     "excluded only by the sequence-matching criteria described in the following sections.")
        text += self.replaced_note("fastqc") + self.replaced_note("multiqc")
        return ["### Read quality assessment", text]

    def libraries_section(self):
        tno = self.table_no()
        test_note = ""
        if self.has_contrasts:
            uses_controls = any(r.ok and (r.opt("--control-sgrna") or r.opt("--control-gene"))
                                for r in self.test_runs.values())
            if not uses_controls:
                test_note = ", not in the statistical analysis"
        one = len(self.libraries) == 1
        text = (f"{'The sgRNA library was' if one else 'sgRNA libraries were'} supplied as "
                f"{'a CSV file' if one else 'CSV files'} (columns gene, gRNA_ID and sequence; Appendix B) and "
                f"converted to MAGeCK's library format ({tno}). sgRNAs whose gene name contains “Control”, “NTC” or "
                f"“Non-targeting” (case-insensitive) were classed as control sgRNAs; this classification was used "
                f"only to show controls separately in quality control{test_note}."
                + self.replaced_note("build_mageck_library"))
        mode = self.info.get("counting_mode")
        rows = []
        for lib in self.libraries:
            seqs = self.lib_seqs.get(lib["name"])
            if seqs is None:
                rows.append([lib["name"], "n/a", "n/a", "n/a", "n/a", ""])
                continue
            controls = [x for x in seqs if CONTROL_PATTERN.search(x[2] or "")]
            genes = {g for _, _, g in seqs if not CONTROL_PATTERN.search(g or "")}
            lengths = self.lib_lengths(lib["name"])
            samples = (f"all {len(lib['samples'])}" if mode == "all_vs_all"
                       else ", ".join(lib["samples"]) if len(lib["samples"]) <= 12
                       else f"{len(lib['samples'])} samples")
            length_text = (f"{lengths[0]}–{lengths[-1]}" if len(lengths) > 1 else
                           str(lengths[0]) if lengths else "n/a")
            rows.append([lib["name"], f"{len(seqs):,}", f"{len(genes):,}", f"{len(controls):,}", length_text, samples])
        blocks = ["### sgRNA libraries", text]
        n = len(self.samples)
        if mode == "all_vs_all":
            blocks.append(
                "No library was assigned to individual samples, so every sample was counted against every "
                f"library ({count_of(len(self.libraries), 'library', 'libraries')}). For each sample, the share of "
                "its reads matching each library, and any library (reads of an sgRNA sequence present in several "
                "libraries counted once), was then tabulated with the pipeline script "
                f"{code('library_breakdown.py')} to identify the library or libraries the sample contained."
                + self.replaced_note("library_breakdown"))
        elif mode == "assigned":
            blocks.append(
                f"Each of the {count_of(n, 'sample')} was counted only against the library assigned to it in the "
                f"samplesheet ({tno}). The share of each sample's reads matching its library was tabulated with "
                f"the pipeline script {code('library_breakdown.py')}." + self.replaced_note("library_breakdown"))
        else:
            blocks.append("The sample was counted against this library." if n == 1 else
                          f"All {n} samples were counted against this library.")
        blocks.append(f"**{tno}.** sgRNA libraries. Genes exclude controls; “Samples” lists the samples "
                      f"counted against each library.\n\n"
                      + table(["Library", "sgRNAs", "Genes", "Control sgRNAs", "sgRNA length (nt)", "Samples"], rows))
        dups = [(lib, r.duplicates) for lib, r in self.count_runs.items() if r.duplicates]
        if dups:
            blocks.append(
                "Some sgRNAs shared their sequence with an earlier sgRNA of the same library ("
                + oxford([f"{k} in {lib}" for lib, k in dups]) + "). MAGeCK keeps only the first sgRNA listed for "
                "each sequence, so all reads of a shared sequence were counted for that sgRNA and the later "
                "duplicates do not appear in the count tables.")
        return blocks

    def positioning(self):
        title = "### Locating the sgRNA in each read"
        if self.position == "fixed":
            return [title] + self._fixed_position()
        return [title] + self._anchor_discovery() + self._anchor_trimming()

    def _fixed_position(self):
        trim5 = str((self.info.get("guide_position") or {}).get("trim5", ""))
        logged = {r.opt("--trim-5") for r in self.count_runs.values() if r.ok}
        if logged:
            trim5 = ",".join(sorted(logged))
        if trim5.upper() == "AUTO":
            autos = {}   # offset(s) -> [(library, FASTQ)]
            for lib, run in self.count_runs.items():
                for fq, val in run.auto_trim.items():
                    autos.setdefault(val, []).append((lib, fq))
            if len(autos) == 1:
                val, pairs = next(iter(autos.items()))
                found = f"{val} nt for every FASTQ file counted ({count_of(len(pairs), 'sample–library combination')})"
            elif autos:
                found = "; ".join(f"{v} nt ({count_of(len(f), 'sample–library combination')})"
                                  for v, f in sorted(autos.items()))
            else:
                found = to_complete("the offsets MAGeCK chose (not found in its log)")
            return [f"The position of the sgRNA within the reads was determined by MAGeCK count for each FASTQ "
                    f"file ({code('--trim-5 AUTO')}): it tested every offset in the first 100,000 reads, recorded "
                    f"for each read the first offset at which the read contained an exact library sgRNA, and used "
                    f"the most frequent offset(s). Offsets used: {found}." + self.replaced_note("mageck_count")]
        return [f"The sgRNA was assumed to start at the same position in every read: MAGeCK count skipped the "
                f"first {trim5} bases of each read ({code('--trim-5 ' + trim5)}), the length of the constant "
                f"sequence between the start of the "
                f"{'demultiplexed read (after barcode removal)' if self.demux else 'read'} and the sgRNA."]

    def _anchor_discovery(self):
        tno = self.table_no()
        d = dict(self.discover_defaults)
        d.update(self._discover_overrides())
        d = {k: as_number(v) for k, v in d.items()}
        n_reads = d.get("n-reads", 200000)
        rows, unrefined, fewer = [], [], False
        for lib in self.libraries:
            name = lib["name"]
            anchor = self.anchors.get(name) or ""
            disc = self.discovery.get(name) or []
            sampled = None
            if disc:   # count / % of reads gives the number examined (to the rounding of the %)
                c, p = to_float(disc[0].get("count")), to_float(disc[0].get("pct_of_reads"))
                if c and p:
                    sampled = int(round(c * 100 / p))
                    if abs(sampled - n_reads) <= 0.001 * n_reads:
                        sampled = n_reads
                    fewer = fewer or sampled < n_reads
            ref = self.refinement.get(name) or []
            chosen = next((r for r in ref if str(r.get("chosen", "")).startswith("True")), None)
            trimmed = chosen.get("bases_trimmed") if chosen else ""
            if chosen and not (to_float(chosen.get("reads_with_library_guide_after")) or 0):
                unrefined.append(name)
                trimmed = "0 (no library sgRNA found)"
            rows.append([name, code(anchor) if anchor else "n/a", len(anchor) if anchor else "n/a",
                         lib["samples"][0] if lib.get("samples") else "",
                         f"{sampled:,}" if sampled else "n/a", trimmed])
        caption = (f"**{tno}.** Anchor (constant vector sequence immediately 5′ of the sgRNA) used for each "
                   f"library, the sample it was identified from, the number of reads examined, and the number of "
                   f"3′ bases removed by library-aware refinement.\n\n"
                   + table(["Library", "Anchor (5′→3′)", "Length (nt)", "Identified from", "Reads examined",
                            "3′ bases removed"], rows))
        if self.replaced("discover_anchor"):
            return [f"The sgRNA was located in each read from the constant vector sequence immediately 5′ of it "
                    f"(the “anchor”), identified for each library ({tno}) with the project-specific script "
                    f"{code(self.scripts['discover_anchor']['path'])}. "
                    + to_complete("how this script identified the anchor"), caption]
        reads = num(n_reads)
        ratio = to_float(d.get("ratio-cutoff", 0.7)) or 0.7
        text = (
            "Because the sgRNA's position within the read depends on the vector and amplicon design (and differs "
            "between reads when staggered primers are used), the sgRNA was located from the constant vector "
            "sequence immediately 5′ of it (the “anchor”) rather than from a fixed offset. The anchor of each "
            f"library was identified from the reads themselves, without prior knowledge of its sequence or length, "
            f"with the pipeline script {code('discover_anchor.py')}, using the first {reads} reads of one sample "
            f"(the first sample counted against the library{', or all of its reads if it had fewer' if fewer else ''}; "
            f"{tno}). Starting at the first base of the "
            f"demultiplexed read, the most frequent read prefix was extended one base at a time (prefix lengths "
            f"{DISCOVERY_START_K}–{num(d.get('max-length', 40))} nt) until the fraction of reads carrying the most "
            f"frequent prefix fell below {ratio * 100:.0f}% of its value at the previous length; the last prefix "
            f"before this drop was taken as the anchor, provided it was at least {num(d.get('min-length', 8))} nt "
            f"long and shared by at least {num(d.get('min-frequency', 30))}% of the reads examined. A base shared "
            f"by nearly all sgRNAs of a library (such as a 5′ G) cannot be told apart from vector sequence this "
            f"way, so the anchor was then checked against the library: versions shortened by 0–"
            f"{num(d.get('max-refine-trim', 8))} nt at the 3′ end were scored by the number of examined reads in "
            f"which they were immediately followed by an exact library sgRNA, and the version with the most such "
            f"reads (the longer on ties) was used.")
        if unrefined:
            text += (f" For {oxford(unrefined)}, no examined read had a library sgRNA after any version, so the "
                     f"anchor was used as first identified.")
        text += self.extra_note("discover_anchor", tool_logged=False)
        return [text, caption]

    def _discover_overrides(self):
        args = (self.extra_args.get("discover_anchor") or "").split()
        out, i = {}, 0
        while i < len(args):
            tok = args[i]
            if tok.startswith("--") and "=" in tok:
                k, v = tok[2:].split("=", 1)
                out[k] = v
            elif tok.startswith("--") and i + 1 < len(args):
                out[tok[2:]] = args[i + 1]
                i += 1
            i += 1
        return out

    def _anchor_trimming(self):
        runs = {lib: r for lib, r in self.trim_runs.items() if r and r.ok}
        if not runs:
            return ["Each read was then trimmed through its library's anchor"
                    + (f" with the project-specific script {code(self.scripts['guide_anchor_trim']['path'])}"
                       if self.replaced("guide_anchor_trim") else "")
                    + ". " + to_complete("how reads were trimmed through the anchor; no cutadapt log was found")]
        first = next(iter(runs.values()))
        signature = {(r.error_rate, r.no_indels, r.discard_untrimmed, r.other_text(), tuple(r.quality_phrases()))
                     for r in runs.values()}
        lengths = sorted({len(s) for r in runs.values() for _, s, _ in r.adapter_seqs()})
        anchored = all(a for r in runs.values() for _, _, a in r.adapter_seqs())
        flags = show(first.opts, {"-e", "--error-rate", "--errors", "--no-indels"})
        # Minimum overlap: the same rule for every library (e.g. anchor length - 2)?
        overlaps = {}
        for lib, r in runs.items():
            seqs = r.adapter_seqs()
            if seqs and r.overlap is not None:
                overlaps[lib] = len(seqs[0][1]) - r.overlap
        if anchored:
            overlap_text = "in full at the start of the read"
        elif first.overlap is None:
            overlap_text = "over at least 3 nt (cutadapt's default minimum overlap)"
        elif len(set(overlaps.values())) == 1:
            slack = next(iter(overlaps.values()))
            values = sorted({r.overlap for r in runs.values()})
            shown = code(f"--overlap {values[0]}") if len(values) == 1 else code("--overlap")
            overlap_text = ("over its full length" if slack <= 0 else
                            f"over at least its length minus {slack} nt ({shown}), so that an anchor "
                            f"truncated by up to {slack} nt at the start of the read was also found")
        else:
            overlap_text = ("over a minimum length of " +
                            oxford([f"{r.overlap} nt ({lib})" for lib, r in runs.items()]))
        text = (f"Each read was then trimmed through its library's anchor with cutadapt"
                f"{self.v('cutadapt', 'GUIDE_ANCHOR_TRIM', first.version)} ({self.cite('martin2011')}): the anchor "
                f"was given as a 5′ adapter ({code('-g <anchor>')}){', anchored at the read start' if anchored else ''} "
                f"and found {overlap_text}, {error_phrase(first.error_rate, first.no_indels, lengths, flags)}. "
                f"The anchor and every base 5′ of it were removed, so that the sgRNA started at the first base of "
                f"each trimmed read, whatever the length of the sequence preceding the anchor. ")
        if first.discard_untrimmed:
            text += (f"Reads in which the anchor was not found were discarded ({code('--discard-untrimmed')}) and "
                     f"so were not part of the mapping-rate denominator.")
        else:
            text += ("Reads in which the anchor was not found were kept unmodified rather than discarded, so that "
                     "every read of the sample remained in the mapping-rate denominator; as their first bases are "
                     "not an sgRNA, these reads rarely matched one and were mostly left uncounted.")
        if first.quality:
            text += f" In addition, cutadapt applied {oxford(first.quality_phrases())}."
        if first.other_text():
            text += f" Other cutadapt options: {code(first.other_text())}."
        if len(signature) > 1:
            text += (" Libraries were not all trimmed with the same options; Appendix A lists the command used "
                     "for each.")
        missing = [lib for lib in self.lib_names if lib not in runs]
        if missing:
            text += f" No trimming log was found for {oxford(missing)}."
        return [text + self.replaced_note("guide_anchor_trim")]

    def counting(self):
        runs = {lib: r for lib, r in self.count_runs.items() if r.ok}
        title = "### sgRNA read counting"
        if not runs:
            return [title, "Reads were counted per library with MAGeCK count" + self.v("MAGeCK", "MAGECK_COUNT")
                    + f" ({self.cite('li2014')})." + self.replaced_note("mageck_count") + " "
                    + to_complete("the counting parameters; no MAGeCK log was found")]
        first = next(iter(runs.values()))
        version = first.version
        lengths = sorted({L for lib in runs for L in self.lib_lengths(lib)})
        trim = first.opt("--trim-5") or "AUTO"
        if self.position == "anchor":
            start = f"from the first base of the anchor-trimmed read ({code('--trim-5 ' + trim)})"
        elif trim.upper() == "AUTO":
            start = "at the offset determined above"
        else:
            start = f"after the first {trim} bases ({code('--trim-5 ' + trim)})"
        len_text = (f"the length of the library's sgRNAs ({oxford(lengths)} nt)" if len(lengths) <= 1 else
                    f"each sgRNA length present in the library ({oxford(lengths)} nt; tried in turn)")
        text = (f"Reads were counted{' for each library separately' if len(self.libraries) > 1 else ''} with "
                f"MAGeCK count{self.v('MAGeCK', 'MAGECK_COUNT', version)} "
                f"({self.cite('li2014')}). For each read, MAGeCK takes the segment that starts {start} and has "
                f"{len_text}, and counts the read for an sgRNA only if this segment is identical to the sgRNA "
                f"sequence: MAGeCK count allows no mismatches.")
        if present(first.opts, "--unmapped-to-file"):
            text += (f" Reads that matched no sgRNA were written to a separate file ({code('--unmapped-to-file')}) "
                     f"and used only for quality control.")
        norm = first.opt("--norm-method") or "median"
        if norm == "median":
            text += (f" MAGeCK also wrote median-normalised counts ({code('--norm-method median')}): each sample's "
                     f"counts were divided by a size factor equal to the median, over sgRNAs with at least one read "
                     f"in the library's samples, of the ratio of the sample's count to the sgRNA's geometric mean "
                     f"count across those samples (computed from counts + 1).")
        else:
            text += f" Normalised counts were computed with {code('--norm-method ' + norm)}."
        text += (" Raw counts were used for quality control"
                 + (" and for statistical testing, which normalises them itself." if self.has_contrasts else "."))
        blocks = [title, text]
        fallback = []
        for lib, r in runs.items():
            if r.total_norm:
                labels, why = r.labels(), {}
                for i, reason in r.total_norm:   # a zero median is the decisive reason when both apply
                    name = labels[i] if i < len(labels) else f"sample {i + 1}"
                    if name not in why or reason.startswith("a median"):
                        why[name] = reason
                fallback.append(f"{lib} ({'; '.join(f'{s}: {w}' for s, w in why.items())})")
        if fallback and norm == "median":
            blocks.append(
                "For some libraries MAGeCK replaced median normalisation by total-count normalisation (every "
                "sample scaled to the mean number of sgRNA-matching reads across the library's samples), because "
                "at least one sample had a median count of zero or more than 45% of sgRNAs with zero reads, which "
                "makes the median ratio unstable. MAGeCK applies this to all samples of the table. Libraries "
                "and the samples that triggered it: " + ", ".join(fallback)
                + ". This affects only the normalised count tables"
                + (" and the normalisation within MAGeCK test, described below." if self.has_contrasts else "."))
        dropped = {lib: self.dropped_samples(lib) for lib in self.lib_names}
        dropped = {k: v for k, v in dropped.items() if v}
        if dropped:
            blocks.append(
                "MAGeCK cannot normalise a sample without any matching read, and in this run such a sample is "
                "expected (every sample was counted against every library). Samples with no read matching a "
                "library were therefore left out of that library's count table, which was counted again without "
                "them: " + "; ".join(f"{lib}: {', '.join(s)}" for lib, s in dropped.items()) + ".")
        for status in self.other_drop_statuses():
            by_lib = {lib: self.dropped_samples(lib, status) for lib in self.lib_names}
            by_lib = {k: v for k, v in by_lib.items() if v}
            blocks.append(
                f"The counting script also left samples out of count tables with the status {code(status)}: "
                + "; ".join(f"{lib}: {', '.join(s)}" for lib, s in by_lib.items()) + ". "
                + to_complete(f"why samples were left out with the status {status}"))
        others = sorted({f"{n} {v}" if v is not None else n for r in runs.values()
                         for n, v in r.other(MAGECK_COUNT_DESCRIBED) if not isinstance(v, list)})
        if others:
            blocks.append(f"Other MAGeCK count options: {code(' '.join(others))}.")
        blocks[1] += self.replaced_note("mageck_count")
        return blocks

    def library_qc(self):
        title = "### Library representation quality control"
        r_ver = self.version("R", "GUIDE_REPRESENTATION_WRANGLE") or ""
        year = re.search(r"\((\d{4})-", r_ver)   # r_ver: "4.3.3 (2024-02-29)"
        r_num = r_ver.split(" ")[0]
        pkgs = []
        for p in R_PACKAGES:
            pv = self.version(f"R package {p}", "GUIDE_REPRESENTATION_WRANGLE")
            pkgs.append(f"{p} v{pv}" if pv else p)
        r_cite = self.cite("rcore", year=year.group(1) if year else "n.d.")
        software = (f"R{' v' + r_num if r_num else ''} ({r_cite}) with the tidyverse packages {oxford(pkgs)} "
                    f"({self.cite('wickham2019')})")
        if self.replaced("guide_representation_wrangle"):
            return [title, f"Library representation was assessed with the project-specific script "
                           f"{code(self.scripts['guide_representation_wrangle']['path'])}, run in {software}. "
                    + to_complete("the metrics this script computes")]
        uneven = [lib for lib in self.lib_names if len(self.lib_lengths(lib)) > 1]
        classified = [lib for lib in self.lib_names if len(self.lib_lengths(lib)) == 1]
        multi = [lib for lib in self.lib_names if len(self.active_samples(lib)) >= 2]
        text = (f"Library representation was assessed for every sample from its raw counts with the pipeline's R "
                f"scripts, run in {software}. Samples without any counted read in a library were excluded from that "
                f"library's assessment. For each sample, the following were computed:")
        items = [
            "**Mapping rate**: the percentage of the sample's reads (as given to MAGeCK) that were counted, i.e. "
            "matched a library sgRNA exactly.",
            f"**Gini index** of sgRNA read counts ({self.cite('li2015')}), computed over all sgRNAs of the library "
            "including those with zero reads (so dropout increases it); 0 indicates perfectly even representation.",
            "**90/10 ratio**: the 90th divided by the 10th percentile of read counts across all library sgRNAs, "
            "including zero-count ones (the 10th percentile floored at 1); flagged when ≥ 10.",
            "**Undetected sgRNAs**: the percentage of library sgRNAs with no reads; flagged when ≥ 5%.",
            "**Representation within 2 log2 units**: the percentage of sgRNAs whose log2(count + 1) lies within "
            "2 of the log2 of the sample's mean count.",
            "**Sequencing depth**: total reads per sample (from the MAGeCK count log), reads matching a library "
            "sgRNA, and reads per sgRNA (reads matched to the library divided by the number of sgRNAs in the library).",
        ]
        if classified:
            note = (f" Not computed for {oxford(uneven)}, whose sgRNAs differ in length." if uneven else "")
            items.append("**Unmapped-read classification**: read segments MAGeCK could not assign "
                         f"({code('--unmapped-to-file')}) that differ from a library sgRNA by exactly one "
                         "substitution (Hamming distance 1) were classed as near-misses (likely sequencing or PCR "
                         "errors) and all others as far-misses." + note)
        if multi:
            items.append("**Concordance between samples**: Pearson correlation of log10(count + 1) and Spearman "
                         "correlation of raw counts, for every pair of samples counted against the same library; the "
                         "samples were ordered in the heatmaps by hierarchical clustering (average linkage on 1 minus "
                         "the Pearson correlation).")
        pca_libs = [lib for lib in self.lib_names if len(self.active_samples(lib)) >= 3]
        if pca_libs:
            note = (" Libraries with fewer than three samples were not analysed."
                    if len(pca_libs) < len(self.lib_names) else "")
            items.append("**Principal component analysis**: sample counts were converted to log2(counts per million "
                         "+ 1) (counts per million = count divided by the sample's total counts, times 10^6), guides "
                         "without reads in any sample were removed, and a principal component analysis was run with "
                         f"R's {code('prcomp')} (centred, not scaled) on the samples; the percentage of variance "
                         "explained by each component is reported." + note)
        plots = (f"Lorenz curves, count distributions (ECDF, with the mean read count per sgRNA of each sample), "
                 f"GC content against read depth, "
                 + ("principal component plots, " if pca_libs else "")
                 + ("correlation heatmaps, " if multi else "")
                 + f"and read depth were plotted with ggplot2 "
                 f"({self.cite('wickham2016')}) and assembled with the metrics into an HTML report for each library"
                 + (", and one report comparing all libraries" if len(self.libraries) > 1 else "") + ".")
        blocks = [title, text, "\n".join(f"- {i}" for i in items), plots]
        for s in ("guide_representation_plot", "guide_representation_report", "guide_representation_compare_report",
                  "my_theme"):
            if self.replaced(s):
                blocks[-1] += self.replaced_note(s)
        return blocks

    def statistics(self):
        if not self.has_contrasts:
            return []
        title = "### Statistical analysis"
        tno = self.table_no()
        runs = {name: r for name, r in self.test_runs.items() if r.ok}
        if not self.contrasts:
            return [title, "Comparisons were requested (contrasts file), but no resolved contrast was found. "
                    + to_complete("the statistical analysis")]
        rows = []
        for c in self.contrasts:
            lib = self._library_of_table(c.get("count_table", ""))
            rows.append([c["name"], ", ".join(c["treat_cols"].split(",")), ", ".join(c["control_cols"].split(",")), lib])
        caption = (f"**{tno}.** Comparisons tested. Each used the count table of the library holding both "
                   f"conditions.\n\n" + table(["Comparison", "Treatment samples", "Control samples", "Library"], rows))
        if not runs:
            return [title, f"Each comparison ({tno}) was tested with MAGeCK test{self.v('MAGeCK', 'MAGECK_TEST')} "
                           f"({self.cite('li2014')})." + self.replaced_note("mageck_test") + " "
                    + to_complete("the test parameters; no MAGeCK test log was found"), caption]
        first = next(iter(runs.values()))

        def setting(name, default):
            vals = {r.opt(name) or default for r in runs.values()}
            return vals.pop() if len(vals) == 1 else None

        norm = setting("--norm-method", "median")
        adjust = setting("--adjust-method", "fdr")
        remove_zero = setting("--remove-zero", "both")
        rz_threshold = setting("--remove-zero-threshold", "0")
        lfc = setting("--gene-lfc-method", "median")
        alpha = setting("--gene-test-fdr-threshold", "0.25")
        paired = {present(r.opts, "--paired") for r in runs.values()}
        varsamples = {r.opt("--variance-estimation-samples") for r in runs.values()}

        def flag(name, value, default):
            if value is None:
                return f"{code(name)} differed between comparisons; see Appendix A"
            return (code(f"{name} {value}") if any(r.opt(name) for r in runs.values())
                    else f"MAGeCK's default, {code(f'{name} {default}')}")

        n_controls = {c["name"]: len([x for x in c["control_cols"].split(",") if x]) for c in self.contrasts}
        if varsamples != {None}:
            var_text = "from the samples given by " + code("--variance-estimation-samples")
        elif all(k >= 2 for k in n_controls.values()):
            var_text = "from the control samples"
        elif all(k < 2 for k in n_controls.values()):
            var_text = "from all samples of the comparison, as each comparison had a single control sample"
        else:
            var_text = ("from the control samples when a comparison had at least two, and otherwise from all of "
                        f"its samples ({tno})")
        norm_text = {"median": "median-normalised as described above",
                     "total": "normalised by total read count",
                     "none": "not normalised",
                     "control": "normalised to the control sgRNAs given"}.get(norm or "", "normalised")
        adjust_text = {"fdr": "the Benjamini–Hochberg procedure", "holm": "Holm's method",
                       "pounds": "Pounds' method"}.get(adjust or "", "the method shown")
        adjust_ref = f"{self.cite('benjamini1995')}; " if adjust == "fdr" else ""
        zero_where = {"both": "in both conditions", "control": "in the control samples",
                      "treatment": "in the treatment samples", "any": "in either condition"}.get(remove_zero or "")
        zero_level = "zero" if rz_threshold in ("0", "0.0") else f"at most {rz_threshold}"
        if remove_zero == "none":
            zero_text = f"No sgRNA was removed for low counts ({code('--remove-zero none')})."
        else:
            zero_text = (f"sgRNAs whose mean normalised count was {zero_level} {zero_where or ''} were removed "
                         f"({flag('--remove-zero', remove_zero, 'both')}).")
        text = (
            f"Genes depleted or enriched in each comparison ({tno}) were identified with MAGeCK test"
            f"{self.v('MAGeCK', 'MAGECK_TEST', first.version)} ({self.cite('li2014')}), using its robust rank "
            f"aggregation (RRA) algorithm on the raw count table of the library holding both conditions. Counts "
            f"were {norm_text} ({flag('--norm-method', norm, 'median')}). MAGeCK estimated the variance of each "
            f"sgRNA's count from its mean, with a mean–variance model fitted across all sgRNAs {var_text}, and "
            f"computed negative-binomial p-values for lower (depletion) and higher (enrichment) counts in treatment "
            f"than in control samples; sgRNA-level p-values were adjusted with {adjust_text} "
            f"({adjust_ref}{flag('--adjust-method', adjust, 'fdr')}). {zero_text} Genes were ranked separately for "
            f"depletion and enrichment with α-RRA: for each gene, only its sgRNAs ranked within the top α fraction of "
            f"all sgRNAs were considered, α being the fraction of sgRNAs with a p-value ≤ {alpha} "
            f"({flag('--gene-test-fdr-threshold', alpha, '0.25')}; bounded to 0.05–0.5 for depletion). Gene-level "
            f"p-values were obtained by permutation and converted into false discovery rates by MAGeCK. Gene log2 "
            f"fold changes are the {lfc} of the gene's sgRNA log2 fold changes "
            f"({flag('--gene-lfc-method', lfc, 'median')}).")
        if any(len(c["treat_cols"].split(",")) > 1 or len(c["control_cols"].split(",")) > 1 for c in self.contrasts):
            text += (" Treatment and control samples were paired in the order listed (" + code("--paired") + ")."
                     if paired == {True} else " Samples of the same condition were treated as unpaired replicates.")
        ctrl = {r.opt("--control-sgrna") or r.opt("--control-gene") for r in runs.values()} - {None}
        if ctrl:
            text += (f" Control sgRNAs ({code(' '.join(sorted(ctrl)))}) were used for normalisation and as the "
                     f"null distribution of RRA.")
        else:
            text += " All sgRNAs, not only control sgRNAs, were used for normalisation and as the RRA null distribution."
        if None in (norm, adjust, remove_zero, rz_threshold, lfc, alpha) or len(paired) > 1:
            text += " Not every comparison was run with the same options; Appendix A lists each command."
        others = sorted({f"{n} {v}" if v is not None else n for r in runs.values()
                         for n, v in r.other(MAGECK_TEST_DESCRIBED) if not isinstance(v, list)})
        if others:
            text += f" Other MAGeCK test options: {code(' '.join(others))}."
        fallback = [name for name, r in runs.items() if r.total_norm]
        if fallback and norm == "median":
            text += (f" For {oxford(fallback)}, MAGeCK replaced median normalisation by total-count normalisation "
                     f"because a sample had a median count of zero or more than 45% of sgRNAs with zero reads.")
        return [title, text + self.replaced_note("mageck_test") + self.replaced_note("validate_contrasts"), caption]

    def _library_of_table(self, table_name):
        for lib in self.lib_names:
            if table_name == f"{self.project}_{lib}.count.txt":
                return lib
        return table_name

    def software(self):
        pipe = self.info.get("pipeline") or {}
        nf = self.info.get("nextflow") or {}
        tno = self.table_no()
        engine = nf.get("container_engine")
        engine_v = re.search(r"(\d+\.\d+[\w.\-]*)", nf.get("container_engine_version") or "")
        commit = pipe.get("commit")
        version = (f"version {pipe.get('version')}" if pipe.get("version") else "") + \
            (f", git commit {commit}" if commit else "") + (" with uncommitted changes" if pipe.get("modified") else "")
        cluster = " on an LSF compute cluster" if "pmacs" in (nf.get("profile") or "") else ""
        text = (f"The analysis was run with the {pipe.get('name') or 'CRISPR_SCREEN'} Nextflow pipeline (Richard "
                f"Phillips Lab; {version}) under Nextflow v{nf.get('version', '?')} ({self.cite('ditommaso2017')}){cluster}. ")
        if engine:
            text += (f"Every step ran in a fixed container image ({tno}) with "
                     f"{'Apptainer' if engine == 'apptainer' else engine.capitalize()}"
                     f"{' v' + engine_v.group(1) if engine_v else ''} ({self.cite('kurtzer2017')}); cutadapt, "
                     f"MAGeCK, FastQC and MultiQC images came from BioContainers ({self.cite('daveiga2017')}). ")
        text += (f"The run's parameters are in {code(nf.get('config_file') or 'config.yaml')}, its resolved "
                 f"configuration in {code(self.params.get('metadata_dir', './metadata').rstrip('/') + '/pipeline_info/run_info.json')}, "
                 f"and the exact command of every step in a {code('*_command.sh')} file next to that step's output.")
        rows = []
        logged_cut = next((r.version for r in [self.demux_run] + list(self.trim_runs.values()) if r and r.version), None)
        logged_mageck = next((r.version for r in self.count_runs.values() if r.version), None)
        for tool, uses in STEP_TOOLS.items():
            ran = [(p, what) for p, what in uses if p in self.processes]
            if not ran:
                continue
            logged = logged_cut if tool == "cutadapt" else logged_mageck if tool == "MAGeCK" else None
            ver = logged or self.version(tool, ran[0][0]) or "not recorded"
            purposes = []
            for _, what in ran:
                if what not in purposes:
                    purposes.append(what)
            rows.append([tool, ver, oxford(purposes), code(self.containers.get(ran[0][0], "n/a"))])
        py = next((p for p in PYTHON_PROCESSES if p in self.processes), None)
        if py:
            rows.append(["Python", self.version("Python", py) or "not recorded",
                         "pipeline helper scripts (library conversion"
                         + (", anchor discovery" if self.position == "anchor" else "")
                         + (", reads per library" if len(self.libraries) > 1 else "")
                         + (", contrast resolution" if self.has_contrasts else "") + ")",
                         code(self.containers.get(py, "n/a"))])
        r_ver = self.version("R", "GUIDE_REPRESENTATION_WRANGLE")
        rows.append(["R", r_ver.split(" ")[0] if r_ver else "not recorded", "library representation QC, plots and reports",
                     code(self.containers.get("GUIDE_REPRESENTATION_WRANGLE", "n/a"))])
        for p in R_PACKAGES:
            rows.append([f"R package {p}", self.version(f"R package {p}", "GUIDE_REPRESENTATION_WRANGLE") or "not recorded",
                         "", ""])
        rows.append(["Nextflow", nf.get("version", "?"), "workflow management", ""])
        if engine:
            rows.append(["Apptainer" if engine == "apptainer" else engine, engine_v.group(1) if engine_v else "?",
                         "container runtime", ""])
        blocks = ["### Software and reproducibility", text,
                  f"**{tno}.** Software used in this run, with the container image each ran in.\n\n"
                  + table(["Software", "Version", "Used for", "Container image"], rows)]
        blocks += self._modifications()
        return blocks

    def _modifications(self):
        scripts_run = []
        for p in self.processes:
            for s in PROCESS_SCRIPTS.get(p, []):
                if s not in scripts_run:
                    scripts_run.append(s)
        rows = []
        for s in scripts_run:
            meta = self.scripts.get(s) or {}
            if meta.get("replaced"):
                rows.append([SCRIPT_LABELS.get(s, s), code(s), f"replaced by {code(meta.get('path'))} "
                             f"(MD5 {meta.get('md5', '?')})"])
            if self.extra_args.get(s):
                rows.append([SCRIPT_LABELS.get(s, s), code(s), f"extra arguments {code(self.extra_args[s])}"])
        if not rows:
            return ["All steps used the pipeline's default scripts, with no extra arguments."]
        tno = self.table_no()
        return [f"This run differed from the pipeline's defaults in the steps listed in {tno}; the text above "
                f"describes what these steps actually ran where the tools logged it.",
                f"**{tno}.** Project-specific changes to the pipeline's default steps.\n\n"
                + table(["Step", "Script", "Change"], rows)]

    def references(self):
        if not self.cited:
            return []
        refs = sorted((REFERENCES[k][1].format(**self.ref_fmt) for k in self.cited), key=str.lower)
        return ["## References", "\n".join(f"- {r}" for r in refs)]

    def appendix_commands(self):
        blocks = ["## Appendix A. Exact commands",
                  "As logged by the tools themselves; long sample lists are abbreviated in angle brackets."]
        if self.demux_run and self.demux_run.ok:
            blocks.append("**Demultiplexing (cutadapt)**\n\n```\n" + self._short_cutadapt(self.demux_run) + "\n```")
        trims = {}   # command -> libraries (libraries with the same anchor share one command)
        for lib, r in self.trim_runs.items():
            if r and r.ok:
                trims.setdefault(self._short_cutadapt(r, per_sample=True), []).append(lib)
        if trims:
            blocks.append("**Anchor trimming (cutadapt; one command per sample, shown for one sample)**\n\n```\n"
                          + "\n".join(f"# {', '.join(libs)}\n{cmd}" for cmd, libs in trims.items()) + "\n```")
        counts = [(lib, r) for lib, r in self.count_runs.items() if r.ok]
        if counts:
            blocks.append("**Read counting (MAGeCK count)**\n\n```\n"
                          + "\n".join(f"# {lib}\n{r.short_command()}" for lib, r in counts) + "\n```")
        tests = [(name, r) for name, r in self.test_runs.items() if r.ok]
        if tests:
            blocks.append("**Statistical testing (MAGeCK test)**\n\n```\n"
                          + "\n".join(f"# {name}\n{r.short_command(max_items=99)}" for name, r in tests) + "\n```")
        if len(blocks) == 2:
            blocks.append("No tool command was recorded for this run.")
        return blocks

    @staticmethod
    def _short_cutadapt(run, per_sample=False) -> str:
        """cutadapt + its arguments; beyond the first two adapters (one per
        sample when demultiplexing), the rest are summarized."""
        argv = list(run.argv)
        if per_sample:   # <label>_trimmed.fastq.gz / <label>.fastq.gz -> <sample>
            out = next((t for t in argv if t.endswith("_trimmed.fastq.gz")), None)
            if out:
                label = os.path.basename(out)[: -len("_trimmed.fastq.gz")]
                argv = [re.sub(rf"(^|/){re.escape(label)}(?=[_.])", r"\1<sample>", t) for t in argv]
        adapters, rest, i = [], [], 0
        while i < len(argv):
            if argv[i] in CUTADAPT_ADAPTER_OPTS and i + 1 < len(argv):
                adapters.append(argv[i:i + 2])
                i += 2
            else:
                rest.append(argv[i])
                i += 1
        parts = ["cutadapt"] + [t for a in adapters[:2] for t in a]
        if len(adapters) > 2:
            parts.append(f"<{len(adapters) - 2} more: -g <sample>=<barcode>>")
        return " ".join(parts + rest)

    def appendix_inputs(self):
        rows = []
        nf = self.info.get("nextflow") or {}
        rows.append(["Configuration", code(nf.get("config_file") or "config.yaml"), ""])
        sheet = self.info.get("samplesheet") or {}
        if sheet:
            rows.append(["Samplesheet", code(sheet.get("path")), sheet.get("md5", "")])
        for lib in self.libraries:
            rows.append([f"Library {lib['name']}", code(lib.get("path")), lib.get("md5", "")])
        contrasts = self.info.get("contrasts") or {}
        if contrasts:
            rows.append(["Contrasts", code(contrasts.get("path")), contrasts.get("md5", "")])
        if self.demux:
            rows.append(["Reads (multiplexed R1)", code(self.input.get("fastq_r1")), ""])
            if self.paired:
                rows.append(["Reads (multiplexed R2)", code(self.input.get("fastq_r2")), ""])
        else:
            rows.append(["Reads (per-sample FASTQs)", code(self.input.get("fastq_dir")), ""])
        return ["## Appendix B. Input files",
                "Paths are as recorded at run time (relative to the project directory for the project's own tables; the raw FASTQ files are listed by the full path where they are stored). MD5 checksums identify the exact input tables used.",
                table(["Input", "Path", "MD5"], rows)]

    # ---- assembly -----------------------------------------------------------------
    def render(self) -> str:
        blocks = []
        for name, fn in [("header", self.header), ("demultiplexing", self.demultiplexing),
                         ("read quality", self.read_quality), ("libraries", self.libraries_section),
                         ("sgRNA position", self.positioning), ("counting", self.counting),
                         ("library QC", self.library_qc), ("statistics", self.statistics),
                         ("software", self.software), ("references", self.references),
                         ("appendix A", self.appendix_commands), ("appendix B", self.appendix_inputs)]:
            try:
                blocks += fn()
            except Exception:   # a documentation bug must not fail the run at its last step
                print(f"ERROR while writing the {name} section:", file=sys.stderr)
                traceback.print_exc()
                blocks.append(f"**[The {name} section could not be generated automatically; see "
                              f"methods_command.sh / the METHODS task's .command.err.]**")
        return "\n\n".join(blocks).rstrip() + "\n"


def discover_defaults(path):
    """argparse defaults of the discover_anchor script this run used: {'n-reads': '200000', ...}."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r"""add_argument\(\s*["']--([\w-]+)["'].*?default=([-0-9_.]+)""", line)
            if m:
                out[m.group(1)] = m.group(2).replace("_", "")
    return out


def as_number(value):
    """'200000' -> 200000, '0.70' -> 0.7, '30.0' -> 30 (for display); other values unchanged."""
    f = to_float(value)
    if f is None:
        return value
    return int(f) if f.is_integer() else f


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-info", required=True)
    ap.add_argument("--evidence-dir", required=True)
    ap.add_argument("--discover-script", default=None)
    ap.add_argument("--output", default="METHODS.md")
    ap.add_argument("--run-info-out", default=None)
    ap.add_argument("--versions-out", default=None)
    args = ap.parse_args()

    with open(args.run_info, encoding="utf-8") as fh:
        info = json.load(fh)
    doc = Methods(info, args.evidence_dir, discover_defaults(args.discover_script))
    text = doc.render()
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"  Saved: {args.output}", file=sys.stderr)

    if args.run_info_out:
        with open(args.run_info_out, "w", encoding="utf-8") as fh:
            json.dump(info, fh, indent=2)
            fh.write("\n")
    if args.versions_out:
        with open(args.versions_out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["image", "tool", "version"])
            for image in sorted(doc.versions):
                for tool, ver in sorted(doc.versions[image].items()):
                    w.writerow([image, tool, ver])


if __name__ == "__main__":
    main()
