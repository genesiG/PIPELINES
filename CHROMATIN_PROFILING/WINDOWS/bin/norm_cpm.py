#!/usr/bin/env python3
"""norm_cpm.py -- CPM "scale factors" (library size only) for ONE normalisation run. Standard library only.

CLI contract (called by the NORM_CPM module; --key value):

  norm_cpm.py --ids id1,id2,...   sample ids of the run (comma and/or space separated; used verbatim, order kept)
              --run <name>        run label (only checked to be non-empty; the caller names the output)
              --out <file>        output table, by convention <run>.cpm.scale_factors.tsv

Output <file>: tab-separated, header sample_id, method, raw_factor, bigwig_multiplier, edger_norm_factor, qc_pass, note.
One row per id: method cpm, raw_factor 1, bigwig_multiplier 1, edger_norm_factor 1, qc_pass TRUE, and the note
  "library size only; bigWig denominator = final BAM mapped reads (deepTools); edgeR lib.size = csaw total".
Exit codes: 0 ok; 2 usage error (no ids, duplicated ids, empty run label).

Algebra: CPM does not rescale beyond the library size. bamCoverage --normalizeUsing RPKM|CPM --scaleFactor 1 divides by the
final-BAM mapped reads; edgeR uses lib.size = csaw total with norm.factors = 1. All factors are 1 (the MA reference line of
this method is exactly 0). Invariant: bigwig_multiplier x edger_norm_factor = 1 for every library.
"""

import argparse
import sys

COLUMNS = ["sample_id", "method", "raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "note"]
NOTE = "library size only; bigWig denominator = final BAM mapped reads (deepTools); edgeR lib.size = csaw total"


def die(msg, status=2):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(status)


def main():
    ap = argparse.ArgumentParser(description="CPM (library size only) scale-factor table.")
    ap.add_argument("--ids", nargs="+", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ids = []
    for v in a.ids:
        ids.extend(x.strip() for x in v.split(",") if x.strip())
    if not ids:
        die("--ids is empty")
    if len(set(ids)) != len(ids):
        die("duplicated sample ids: %s" % ", ".join(sorted({i for i in ids if ids.count(i) > 1})))
    if not a.run.strip():
        die("--run is empty")
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sid in ids:
            fh.write("\t".join([sid, "cpm", "1", "1", "1", "TRUE", NOTE]) + "\n")
    print("Wrote %s (%d libraries, run %s)" % (a.out, len(ids), a.run))


if __name__ == "__main__":
    main()
