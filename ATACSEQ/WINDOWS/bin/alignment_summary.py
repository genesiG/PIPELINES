#!/usr/bin/env python3
"""alignment_summary.py -- bowtie2 alignment rates per sample, read from the bowtie2 logs of ALIGN.

Called by ALIGNMENT_SUMMARY (modules/local/alignment_summary.nf):
    alignment_summary.py --logs-dir . --min-pct 20 --out alignment_summary.tsv

Arguments:
  --logs-dir  directory with the <S>.bowtie2.log files (bowtie2 stderr of ALIGN); default "."
  --min-pct   a sample whose overall_alignment_rate is below this percentage is flagged; default 20
  --out       output table; default alignment_summary.tsv

Parsing is ported from the legacy step_2b_alignment_summary.py (paired-end "aligned concordantly" lines; the
single-end "aligned 0 times" wording is the fallback). Values are written as bowtie2 printed them.
Output: alignment_summary.tsv (tab separated), one row per sample, sorted by sample_id:
  sample_id, total_reads, aligned_0_times, aligned_0_times_pct, aligned_1_time, aligned_1_time_pct,
  aligned_multiple_times, aligned_multiple_times_pct, overall_alignment_rate, below_min_pct
  sample_id = file name minus ".bowtie2.log". below_min_pct is "true" when overall_alignment_rate < --min-pct,
  "false" when not, and "NA" when the overall rate could not be parsed (a WARNING is printed for that sample).
  Any other value that cannot be parsed is written as NA.
Warnings: one line per flagged sample on stderr. A low rate never fails the run (exit 0).
Exit 3: no *.bowtie2.log file in --logs-dir.
"""
import argparse
import os
import re
import sys

SUFFIX = ".bowtie2.log"
COLS = ["sample_id", "total_reads", "aligned_0_times", "aligned_0_times_pct", "aligned_1_time",
        "aligned_1_time_pct", "aligned_multiple_times", "aligned_multiple_times_pct",
        "overall_alignment_rate", "below_min_pct"]

TOTAL = r"(\d+) reads; of these:"
RATE = r"([\d.]+)% overall alignment rate"
# (column, pattern); each pattern has the count in group 1 and the percentage in group 2
PE = {
    "aligned_0_times": r"(\d+) \(([\d.]+)%\) aligned concordantly 0 times",
    "aligned_1_time": r"(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time",
    "aligned_multiple_times": r"(\d+) \(([\d.]+)%\) aligned concordantly >1 times",
}
SE = {
    "aligned_0_times": r"(\d+) \(([\d.]+)%\) aligned 0 times",
    "aligned_1_time": r"(\d+) \(([\d.]+)%\) aligned exactly 1 time",
    "aligned_multiple_times": r"(\d+) \(([\d.]+)%\) aligned >1 times",
}


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def parse_log(path, min_pct, sample_id):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    rec = {c: "NA" for c in COLS}
    rec["sample_id"] = sample_id

    m = re.search(TOTAL, text)
    if m:
        rec["total_reads"] = m.group(1)

    pats = PE if re.search(PE["aligned_0_times"], text) else SE
    for col, pat in pats.items():
        m = re.search(pat, text)
        if m:
            rec[col] = m.group(1)
            rec[col + "_pct"] = m.group(2)

    m = re.search(RATE, text)
    if m:
        rec["overall_alignment_rate"] = m.group(1)
        rec["below_min_pct"] = "true" if float(m.group(1)) < min_pct else "false"
    else:
        sys.stderr.write("WARNING: could not read the overall alignment rate of %s\n" % sample_id)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs-dir", default=".")
    ap.add_argument("--min-pct", type=float, default=20.0)
    ap.add_argument("--out", default="alignment_summary.tsv")
    a = ap.parse_args()

    files = sorted(f for f in os.listdir(a.logs_dir) if f.endswith(SUFFIX))
    if not files:
        die("no %s files in %s" % (SUFFIX, a.logs_dir))

    rows = []
    for f in files:
        sid = f[:-len(SUFFIX)]
        rec = parse_log(os.path.join(a.logs_dir, f), a.min_pct, sid)
        if rec["below_min_pct"] == "true":
            sys.stderr.write("WARNING: %s has an overall alignment rate of %s%% (below %g%%)\n"
                             % (sid, rec["overall_alignment_rate"], a.min_pct))
        rows.append(rec)

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLS) + "\n")
        for rec in rows:
            fh.write("\t".join(rec[c] for c in COLS) + "\n")
    print("wrote %s (%d samples)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
