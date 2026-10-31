#!/usr/bin/env python3
"""qc_peak_stats.py -- summary statistics of every peak BED staged in the task directory.

Called by QC_PEAK_STATS (modules/local/qc_peak_stats.nf):
    qc_peak_stats.py --peaks-dir . --out peak_stats.tsv

Arguments:
  --peaks-dir  directory scanned for peak files (default "."): *.bed, *.narrowPeak, *.broadPeak (symlinks followed)
  --out        output table (default peak_stats.tsv)
Output: peak_stats.tsv, one row per peak file: peak_set (file name without extension), file, n_peaks, total_bp,
  min_width, median_width, mean_width, max_width, iqr_width, n_chroms, pct_autosomal (peaks on contigs named
  [chr]<number>; a contig-name test, no genome is assumed), peak_type (narrowPeak | broadPeak | bed, from the file extension).
  Header lines (#, track, browser) are skipped.
  A file without peaks gets n_peaks = 0 and NA for the width columns (not an error).
Exit 3: no peak file found, or a line with fewer than 3 columns / non-integer coordinates.
"""
import argparse
import os
import re
import sys

AUTOSOME = re.compile(r"^(chr)?[0-9]+$", re.IGNORECASE)
EXTS = (".bed", ".narrowPeak", ".broadPeak")
PEAK_TYPE = {".narrowPeak": "narrowPeak", ".broadPeak": "broadPeak"}


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def quantile(sorted_vals, q):
    """Linear-interpolation quantile of a sorted list."""
    if not sorted_vals:
        return float("nan")
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def stats(path):
    widths, chroms = [], {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                die("%s line %d has fewer than 3 tab-separated columns" % (path, ln))
            try:
                w = int(f[2]) - int(f[1])
            except ValueError:
                die("%s line %d has non-integer coordinates" % (path, ln))
            widths.append(w)
            chroms[f[0]] = chroms.get(f[0], 0) + 1
    n = len(widths)
    if n == 0:
        return [0, 0] + ["NA"] * 5 + [0, "NA"]
    widths.sort()
    auto = sum(c for k, c in chroms.items() if AUTOSOME.match(k))
    return [n, sum(widths), min(widths), "%.1f" % quantile(widths, 0.5), "%.1f" % (sum(widths) / n), max(widths),
            "%.1f" % (quantile(widths, 0.75) - quantile(widths, 0.25)), len(chroms), "%.2f" % (100.0 * auto / n)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peaks-dir", default=".")
    ap.add_argument("--out", default="peak_stats.tsv")
    a = ap.parse_args()
    files = sorted(f for f in os.listdir(a.peaks_dir) if f.endswith(EXTS) and os.path.isfile(os.path.join(a.peaks_dir, f)))
    if not files:
        die("no peak files (*.bed, *.narrowPeak, *.broadPeak) in %s" % a.peaks_dir)
    header = ["peak_set", "file", "n_peaks", "total_bp", "min_width", "median_width", "mean_width", "max_width",
              "iqr_width", "n_chroms", "pct_autosomal", "peak_type"]
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for f in files:
            row = stats(os.path.join(a.peaks_dir, f))
            fh.write("\t".join([os.path.splitext(f)[0], f] + [str(x) for x in row] + [PEAK_TYPE.get(os.path.splitext(f)[1], "bed")]) + "\n")
    print("wrote %s (%d peak files)" % (a.out, len(files)))


if __name__ == "__main__":
    main()
