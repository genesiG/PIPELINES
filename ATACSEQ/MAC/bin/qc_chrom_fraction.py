#!/usr/bin/env python3
"""qc_chrom_fraction.py -- per sample, the percentage of raw mapped reads on each chosen contig (e.g. chrM).

Called by QC_CHROM_FRACTION (modules/local/qc_chrom_fraction.nf):
    qc_chrom_fraction.py --input-dir . --chroms=chrM --out chrom_fraction.tsv

Arguments:
  --input-dir  directory with the <S>.raw.idxstats.txt files (BAM_FILTER: samtools idxstats of the unfiltered,
               sorted BAM); default "."
  --chroms     comma-separated contig names (config remove_chroms); empty = no contig columns
  --out        output table; default chrom_fraction.tsv

Input: samtools idxstats format (name, length, mapped, unmapped); the '*' line is ignored.
Output: chrom_fraction.tsv (tab separated), one row per sample, sorted by sample_id:
  sample_id, raw_mapped_reads, then for each requested contig, in the order given: <contig>_reads, <contig>_pct.
  raw_mapped_reads = sum of mapped reads over all contigs (the % denominator).
  <contig>_reads = mapped reads on that contig; <contig>_pct = 100 * <contig>_reads / raw_mapped_reads.
  sample_id = file name minus ".raw.idxstats.txt".
  A requested contig absent from a sample's idxstats: WARNING on stderr, reported as 0.
Exit 3: no *.raw.idxstats.txt in --input-dir, a malformed idxstats line, or a sample with no mapped reads.
"""
import argparse
import os
import sys

SUFFIX = ".raw.idxstats.txt"


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def read_idxstats(path):
    mapped = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                die("%s is not an idxstats file (line: %s)" % (path, line.rstrip("\n")))
            if f[0] == "*":
                continue
            try:
                mapped[f[0]] = int(f[2])
            except ValueError:
                die("%s has a non-integer mapped count (line: %s)" % (path, line.rstrip("\n")))
    return mapped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default=".")
    ap.add_argument("--chroms", default="")
    ap.add_argument("--out", default="chrom_fraction.tsv")
    a = ap.parse_args()

    chroms = [c.strip() for c in a.chroms.split(",") if c.strip()]
    files = sorted(f for f in os.listdir(a.input_dir) if f.endswith(SUFFIX))
    if not files:
        die("no %s files in %s" % (SUFFIX, a.input_dir))

    header = ["sample_id", "raw_mapped_reads"]
    for c in chroms:
        header += [c + "_reads", c + "_pct"]

    rows = []
    for f in files:
        sid = f[:-len(SUFFIX)]
        mapped = read_idxstats(os.path.join(a.input_dir, f))
        total = sum(mapped.values())
        if total <= 0:
            die("%s has no mapped reads; the chromosome fraction is undefined" % f)
        row = [sid, str(total)]
        for c in chroms:
            if c in mapped:
                n = mapped[c]
            else:
                sys.stderr.write("WARNING: contig %s is not in %s; reported as 0\n" % (c, f))
                n = 0
            row += [str(n), "%.4f" % (100.0 * n / total)]
        rows.append(row)

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print("wrote %s (%d samples)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
