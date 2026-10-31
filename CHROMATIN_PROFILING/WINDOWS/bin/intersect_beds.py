#!/usr/bin/env python3
"""intersect_beds.py -- which records of BED A overlap BED B (bedtools intersect -a A -b B -wa -u), or the overlap pieces.

Called by process INTERSECT (intersects.csv row = one call):
    intersect_beds.py --name NAME --bed-a A --bed-b B --label-a LA --label-b LB \
        [--min-overlap 0.05] [--reciprocal true|false] [--output records|overlap]

Arguments:
  --name         output prefix NAME
  --bed-a        BED A; its records are reported
  --bed-b        BED B; used for the overlap test only
  --label-a      label of A (summary)
  --label-b      label of B (summary)
  --min-overlap  minimal overlap fraction of A, 0-1 (bedtools -f F). Empty = no -f (any overlap). Default 0.05
  --reciprocal   true|false (bedtools -r; needs -f, so needs --min-overlap). Default false
  --output       records|overlap: what NAME.intersect.bed holds (see Commands). Default records
Commands (F and -r added only when --min-overlap is set; -r only when --reciprocal true):
  intersect, --output records:  bedtools intersect -a A -b B -wa -u [-f F] [-r]
  intersect, --output overlap:  bedtools intersect -a A -b B [-f F] [-r]   (legacy: one line per A/B overlap, A columns clipped to the overlap)
  only A:     bedtools intersect -a A -b B -v  [-f F] [-r]   (same -f/-r, so the two A sets partition A)
  B overlap:  bedtools intersect -a B -b A -wa -u [-f F] [-r]
Outputs (working directory):
  NAME.intersect.bed          records: A records overlapping B, original columns, each once
                              overlap: overlap pieces, one line per overlapping A/B pair, start/end clipped to the overlap
  NAME.only_a.bed             A records not overlapping B
  NAME.intersect_summary.tsv  header and one row with the columns:
    name label_a label_b n_a n_b n_a_overlapping_b n_a_only pct_a_overlapping_b n_b_overlapping_a min_overlap reciprocal output n_pieces
    n_* are record counts (non-empty lines that are not '#' comments or 'track' lines); the same in both modes (A records).
    n_pieces = lines of NAME.intersect.bed (= n_a_overlapping_b with --output records).
    pct_a_overlapping_b = 100 * n_a_overlapping_b / n_a, 2 decimals.
    min_overlap is NA when --min-overlap is empty.
Exit 3: bed_a or bed_b missing or without records.
Exit 2: --min-overlap not a number in 0-1, --reciprocal not true/false, --reciprocal true without --min-overlap,
        or --output not records/overlap.
Exit 1: bedtools not on PATH, or bedtools failed.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

TRACK_RE = re.compile(r"(track|browser)(\s|$)")


def die(msg, code=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(code)


def count_lines(lines):
    n = 0
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#") or TRACK_RE.match(s):
            continue
        n += 1
    return n


def count_records(path):
    with open(path, encoding="utf-8") as fh:
        return count_lines(fh)


def bedtools_cmd(a, b, mode, min_overlap, reciprocal):
    cmd = ["bedtools", "intersect", "-a", a, "-b", b]
    if mode == "intersect":
        cmd += ["-wa", "-u"]
    elif mode == "pieces":
        pass
    else:
        cmd += ["-v"]
    if min_overlap is not None:
        cmd += ["-f", min_overlap]
        if reciprocal:
            cmd += ["-r"]
    return cmd


def run_to_file(cmd, out_path):
    try:
        with open(out_path, "w", encoding="utf-8") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
    except subprocess.CalledProcessError:
        die("bedtools failed: %s" % " ".join(cmd), 1)


def run_to_text(cmd):
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, check=True, universal_newlines=True)
    except subprocess.CalledProcessError:
        die("bedtools failed: %s" % " ".join(cmd), 1)
    return res.stdout.splitlines()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--bed-a", required=True)
    ap.add_argument("--bed-b", required=True)
    ap.add_argument("--label-a", required=True)
    ap.add_argument("--label-b", required=True)
    ap.add_argument("--min-overlap", default="0.05")
    ap.add_argument("--reciprocal", default="false")
    ap.add_argument("--output", default="records")
    a = ap.parse_args()

    if a.output not in ("records", "overlap"):
        die("--output must be records or overlap (got '%s')" % a.output, 2)

    if a.reciprocal not in ("true", "false"):
        die("--reciprocal must be true or false (got '%s')" % a.reciprocal, 2)
    reciprocal = a.reciprocal == "true"
    mo = a.min_overlap.strip()
    min_overlap = None
    if mo:
        try:
            v = float(mo)
        except ValueError:
            die("--min-overlap must be a fraction 0-1 or empty (got '%s')" % mo, 2)
        if not 0 <= v <= 1:
            die("--min-overlap must be a fraction 0-1 (got '%s')" % mo, 2)
        min_overlap = mo
    if reciprocal and min_overlap is None:
        die("--reciprocal true needs --min-overlap (bedtools -r requires -f)", 2)

    if shutil.which("bedtools") is None:
        die("bedtools not found on PATH", 1)
    for label, path in (("bed_a", a.bed_a), ("bed_b", a.bed_b)):
        if not os.path.isfile(path):
            die("%s not found: %s" % (label, path))
    n_a = count_records(a.bed_a)
    n_b = count_records(a.bed_b)
    if n_a == 0:
        die("bed_a is empty: %s" % a.bed_a)
    if n_b == 0:
        die("bed_b is empty: %s" % a.bed_b)

    out_int = a.name + ".intersect.bed"
    out_only = a.name + ".only_a.bed"
    run_to_file(bedtools_cmd(a.bed_a, a.bed_b, "intersect", min_overlap, reciprocal), out_int)
    run_to_file(bedtools_cmd(a.bed_a, a.bed_b, "only", min_overlap, reciprocal), out_only)
    n_ov = count_records(out_int)
    n_pieces = n_ov
    if a.output == "overlap":
        # counts above stay record-based; NAME.intersect.bed is then replaced by the overlap pieces
        run_to_file(bedtools_cmd(a.bed_a, a.bed_b, "pieces", min_overlap, reciprocal), out_int + ".pieces.tmp")
        os.replace(out_int + ".pieces.tmp", out_int)
        n_pieces = count_records(out_int)
    n_only = count_records(out_only)
    if n_ov + n_only != n_a:
        die("internal: %d overlapping + %d only-A records != %d records in bed_a" % (n_ov, n_only, n_a), 1)
    n_b_ov = count_lines(run_to_text(bedtools_cmd(a.bed_b, a.bed_a, "intersect", min_overlap, reciprocal)))

    pct = "%.2f" % (100.0 * n_ov / n_a)
    header = ["name", "label_a", "label_b", "n_a", "n_b", "n_a_overlapping_b", "n_a_only",
              "pct_a_overlapping_b", "n_b_overlapping_a", "min_overlap", "reciprocal", "output", "n_pieces"]
    row = [a.name, a.label_a, a.label_b, n_a, n_b, n_ov, n_only, pct, n_b_ov,
           mo if mo else "NA", a.reciprocal, a.output, n_pieces]
    with open(a.name + ".intersect_summary.tsv", "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(str(x) for x in row) + "\n")
    print("%s: %d of %d A records overlap B (%s%%)" % (a.name, n_ov, n_a, pct))


if __name__ == "__main__":
    main()
