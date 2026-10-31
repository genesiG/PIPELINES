#!/usr/bin/env python3
"""alignment_summary.py -- HISAT2 paired-end alignment rates per sample, read from the HISAT2 summary files of ALIGN.

Called by ALIGNMENT_SUMMARY (modules/local/alignment_summary.nf):
    alignment_summary.py --logs-dir . --min-pct 50 --out alignment_summary.tsv

Arguments:
  --logs-dir  directory with the <S>.hisat2_summary.txt files (HISAT2 --summary-file of ALIGN); default "."
  --min-pct   a sample whose overall_alignment_rate is below this percentage is flagged; default 50
  --out       output table; default alignment_summary.tsv

Parsing follows the legacy step_2b_alignment_summary.py (GG and AV RNA-seq): the paired-end "aligned concordantly"
lines, the "N reads; of these:" total and the "overall alignment rate" line. The HISAT2 layout is the same as the
bowtie2 one for these lines. Discordant lines ("aligned discordantly 1 time") go only to the aligned_discordantly_1
columns; they are never mixed into the concordant columns. With --no-discordant (ALIGN) that line is normally absent
and its columns are NA. Values are written as HISAT2 printed them.

Output: alignment_summary.tsv (tab separated), one row per sample, sorted by sample_id:
  sample_id, total_pairs, aligned_concordantly_0, aligned_concordantly_0_pct, aligned_concordantly_1,
  aligned_concordantly_1_pct, aligned_concordantly_multi, aligned_concordantly_multi_pct,
  aligned_discordantly_1, aligned_discordantly_1_pct, overall_alignment_rate, below_min_pct
  sample_id = file name minus ".hisat2_summary.txt". below_min_pct is "true" when overall_alignment_rate < --min-pct,
  "false" otherwise.
Warnings: one line per flagged sample on stderr. A low rate never fails the run (exit 0).
Exit 3: no *.hisat2_summary.txt file in --logs-dir; a file without the total, the concordant lines or the
  overall alignment rate line (unparsable summary).
"""
import argparse
import os
import re
import sys

SUFFIX = ".hisat2_summary.txt"
COLS = ["sample_id", "total_pairs", "aligned_concordantly_0", "aligned_concordantly_0_pct",
        "aligned_concordantly_1", "aligned_concordantly_1_pct", "aligned_concordantly_multi",
        "aligned_concordantly_multi_pct", "aligned_discordantly_1", "aligned_discordantly_1_pct",
        "overall_alignment_rate", "below_min_pct"]

TOTAL = r"(\d+) reads; of these:"
RATE = r"([\d.]+)% overall alignment rate"
# (column, pattern); the count is group 1 and the percentage group 2
CONCORDANT = {
    "aligned_concordantly_0": r"(\d+) \(([\d.]+)%\) aligned concordantly 0 times",
    "aligned_concordantly_1": r"(\d+) \(([\d.]+)%\) aligned concordantly exactly 1 time",
    "aligned_concordantly_multi": r"(\d+) \(([\d.]+)%\) aligned concordantly >1 times",
}
DISCORDANT = {
    "aligned_discordantly_1": r"(\d+) \(([\d.]+)%\) aligned discordantly 1 time",
}


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def parse_summary(path, min_pct, sample_id):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    rec = {c: "NA" for c in COLS}
    rec["sample_id"] = sample_id

    m = re.search(TOTAL, text)
    if not m:
        die("%s: no 'reads; of these' line (unparsable HISAT2 summary)" % path)
    rec["total_pairs"] = m.group(1)

    for col, pat in CONCORDANT.items():
        m = re.search(pat, text)
        if not m:
            die("%s: no '%s' line (unparsable HISAT2 summary)" % (path, col))
        rec[col] = m.group(1)
        rec[col + "_pct"] = m.group(2)

    for col, pat in DISCORDANT.items():
        m = re.search(pat, text)
        if m:
            rec[col] = m.group(1)
            rec[col + "_pct"] = m.group(2)

    m = re.search(RATE, text)
    if not m:
        die("%s: no overall alignment rate line (unparsable HISAT2 summary)" % path)
    rec["overall_alignment_rate"] = m.group(1)
    rec["below_min_pct"] = "true" if float(m.group(1)) < min_pct else "false"
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--logs-dir", default=".")
    ap.add_argument("--min-pct", type=float, default=50.0)
    ap.add_argument("--out", default="alignment_summary.tsv")
    a = ap.parse_args()

    files = sorted(f for f in os.listdir(a.logs_dir) if f.endswith(SUFFIX))
    if not files:
        die("no %s files in %s" % (SUFFIX, a.logs_dir))

    rows = []
    for f in files:
        sid = f[:-len(SUFFIX)]
        rec = parse_summary(os.path.join(a.logs_dir, f), a.min_pct, sid)
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
