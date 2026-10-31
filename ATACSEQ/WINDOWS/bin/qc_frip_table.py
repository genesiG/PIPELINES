#!/usr/bin/env python3
"""qc_frip_table.py -- merge the per-sample FRiP files into frip.tsv.

Called by QC_FRIP_TABLE (modules/local/qc_frip_table.nf):
    qc_frip_table.py --input-dir . --samplesheet samplesheet.csv --out frip.tsv

Arguments:
  --input-dir    directory with the <S>.<peak_source>.frip.tsv files from qc_frip.sh (default ".")
  --samplesheet  samplesheet.csv; column sample_id is required. condition, replicate, tmm_group and peak_group are
                 copied from the sheet; a column missing from the sheet is filled with blank (metadata only)
  --out          output table (default frip.tsv)
Output: frip.tsv, one row per (sample, peak set): sample_id, condition, replicate, tmm_group, peak_group, peak_source,
  peaks_file, n_peaks, total_mapped_alignments, reads_in_peaks, frip. Sample and peak source come from the
  CONTENT of each file, never from its name.
Exit 3: no input files, a malformed input file, or a sample_id that is not in the samplesheet.
"""
import argparse
import csv
import os
import sys

META = ["sample_id", "condition", "replicate", "tmm_group", "peak_group"]
FRIP = ["peak_source", "peaks_file", "n_peaks", "total_mapped_alignments", "reads_in_peaks", "frip"]


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default=".")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", default="frip.tsv")
    a = ap.parse_args()

    with open(a.samplesheet, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        ss_rows = list(reader)
    if "sample_id" not in cols:
        die("samplesheet lacks column sample_id")
    ss = {r["sample_id"]: r for r in ss_rows}

    files = sorted(f for f in os.listdir(a.input_dir) if f.endswith(".frip.tsv") and f != os.path.basename(a.out))
    if not files:
        die("no *.frip.tsv in %s" % a.input_dir)
    rows = []
    for f in files:
        with open(os.path.join(a.input_dir, f), newline="", encoding="utf-8") as fh:
            recs = list(csv.DictReader(fh, delimiter="\t"))
        if len(recs) != 1 or any(c not in recs[0] for c in ["sample_id"] + FRIP):
            die("%s is not a one-row FRiP file" % f)
        r = recs[0]
        if r["sample_id"] not in ss:
            die("sample %s (from %s) is not in the samplesheet" % (r["sample_id"], f))
        meta = [r["sample_id"]] + [(ss[r["sample_id"]].get(c) or "") if c in cols else "" for c in META[1:]]
        rows.append(meta + [r[c] for c in FRIP])
    rows.sort(key=lambda x: (x[0], x[5]))
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(META + FRIP) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("wrote %s (%d rows)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
