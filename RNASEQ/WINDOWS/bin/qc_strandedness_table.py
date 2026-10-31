#!/usr/bin/env python3
"""qc_strandedness_table.py -- one table of RSeQC infer_experiment results per sample, with the strandedness call.

Called by QC_STRANDEDNESS_TABLE (modules/local/qc_strandedness_table.nf):
    qc_strandedness_table.py --dir . --strandedness unstranded --threshold 0.8 --out strandedness.tsv

Arguments:
  --dir            directory with the <S>.infer_experiment.txt files (INFER_STRANDEDNESS outputs); default "."
  --strandedness   configured strandedness: unstranded, forward or reverse (config strandedness)
  --threshold      call threshold, in (0.5, 1] (config strand_call_threshold, default 0.8)
  --out            output table; default strandedness.tsv

Per sample (sample_id = file name minus ".infer_experiment.txt"):
  fraction_forward      explained by "1++,1--,2+-,2-+"
  fraction_reverse      explained by "1+-,1-+,2++,2--"
  fraction_undetermined "Fraction of reads failed to determine"
  call                  forward if fraction_forward >= threshold; reverse if fraction_reverse >= threshold;
                        otherwise unstranded
  configured            the --strandedness value
  mismatch              true when call != configured, else false
Fractions are written as RSeQC printed them (4 decimals).

A mismatch prints one "WARN:" line per sample on stderr and never fails the run.
Output: strandedness.tsv (tab separated), one row per sample, sorted by sample_id.
Exit 3: no *.infer_experiment.txt file in --dir, or a file without one of the three fraction lines (unparsable).
"""
import argparse
import glob
import os
import re
import sys

SUFFIX = ".infer_experiment.txt"
COLS = ["sample_id", "fraction_forward", "fraction_reverse", "fraction_undetermined",
        "call", "configured", "mismatch"]
STRANDEDNESS = ("unstranded", "forward", "reverse")
PATTERNS = {
    "undetermined": re.compile(r"Fraction of reads failed to determine:\s*([\d.]+)"),
    "forward": re.compile(r'Fraction of reads explained by "1\+\+,1--,2\+-,2-\+":\s*([\d.]+)'),
    "reverse": re.compile(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s*([\d.]+)'),
}


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def parse_file(path):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    out = {}
    for key, pat in PATTERNS.items():
        m = pat.search(text)
        if not m:
            die("%s: no fraction line for '%s' (unparsable infer_experiment output)" % (path, key))
        out[key] = m.group(1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".")
    ap.add_argument("--strandedness", required=True)
    ap.add_argument("--threshold", type=float, default=0.8)
    ap.add_argument("--out", default="strandedness.tsv")
    a = ap.parse_args()

    if a.strandedness not in STRANDEDNESS:
        die("--strandedness must be one of %s (got '%s')" % (", ".join(STRANDEDNESS), a.strandedness))
    if not (0.5 < a.threshold <= 1.0):
        die("--threshold must be in (0.5, 1] (got %g)" % a.threshold)

    files = sorted(glob.glob(os.path.join(a.dir, "*" + SUFFIX)))
    if not files:
        die("no *%s files in %s" % (SUFFIX, a.dir))

    rows = []
    for path in files:
        sid = os.path.basename(path)[:-len(SUFFIX)]
        fr = parse_file(path)
        fwd, rev = float(fr["forward"]), float(fr["reverse"])
        if fwd >= a.threshold:
            call = "forward"
        elif rev >= a.threshold:
            call = "reverse"
        else:
            call = "unstranded"
        mismatch = call != a.strandedness
        if mismatch:
            sys.stderr.write("WARN: %s: strand call '%s' (forward %s, reverse %s) differs from the configured "
                             "strandedness '%s'\n" % (sid, call, fr["forward"], fr["reverse"], a.strandedness))
        rows.append([sid, fr["forward"], fr["reverse"], fr["undetermined"], call, a.strandedness,
                     "true" if mismatch else "false"])

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLS) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print("wrote %s (%d samples)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
