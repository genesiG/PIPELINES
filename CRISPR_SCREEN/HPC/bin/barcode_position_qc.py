#!/usr/bin/env python3
"""
barcode_position_qc.py
========================
Streams a cutadapt `--info-file` (one line per read; several GB for a full
lane) into a small per-(sample, position) histogram CSV of where in the
read each sample's barcode was actually matched. Plotting is done
separately by bin/plot_barcode_position.R (ggplot2) — this script only
does the memory-bounded data reduction, so it has no plotting dependency
and can run in the same minimal cutadapt container as the demux step
itself.

Info-file columns for a matched 5' adapter: read_id, num_errors,
match_start, match_end, seq_before, seq_matched, seq_after, adapter_name,
qual_before, qual_match, qual_after. Unmatched reads have only 4 columns
(read_id, "-1", sequence, qualities) and are skipped.
"""

import argparse
import csv
import sys
from collections import defaultdict


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--info-file", required=True, help="cutadapt --info-file output")
    ap.add_argument("--project-name", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--max-position", type=int, default=40,
                     help="positions beyond this are clipped into the last bin (default: 40bp)")
    args = ap.parse_args()

    # sample -> {position: count}; kept in memory (small: samples x positions),
    # the info-file itself is streamed line-by-line and never held in memory.
    counts: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    totals: dict[str, int] = defaultdict(int)
    n_matched = 0

    with open(args.info_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue  # unmatched read (discarded, not in any sample's output)
            try:
                start = int(fields[2])
            except ValueError:
                continue
            sample = fields[7]
            if not sample:
                continue
            pos = min(start, args.max_position)
            counts[sample][pos] += 1
            totals[sample] += 1
            n_matched += 1

    if not totals:
        print(f"[WARNING] No matched reads found in {args.info_file} — skipping.", file=sys.stderr)
        sys.exit(0)

    csv_path = f"{args.outdir}/{args.project_name}_barcode_position_counts.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "position", "count", "density"])
        for sample in sorted(totals):
            total = totals[sample]
            for pos in range(args.max_position + 1):
                c = counts[sample].get(pos, 0)
                writer.writerow([sample, pos, c, round(c / total, 5) if total else 0])

    print(f"Saved: {csv_path}")
    print(f"Parsed {n_matched:,} matched reads across {len(totals)} sample(s).")
    for sample in sorted(totals):
        peak_pos = max(range(args.max_position + 1), key=lambda p: counts[sample].get(p, 0))
        peak_pct = round(counts[sample].get(peak_pos, 0) / totals[sample] * 100, 1)
        print(f"  {sample}: peak @ {peak_pos}bp ({peak_pct}% of reads)")


if __name__ == "__main__":
    main()
