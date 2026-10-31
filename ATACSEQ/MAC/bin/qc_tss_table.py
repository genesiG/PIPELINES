#!/usr/bin/env python3
"""qc_tss_table.py -- merge the per-sample TSS enrichment files into tss_enrichment.tsv.

Called by QC_TSS_TABLE (modules/local/qc_tss_table.nf):
    qc_tss_table.py --input-dir . --samplesheet samplesheet.csv --out tss_enrichment.tsv

Arguments:
  --input-dir    directory with the <S>.tss_enrichment.tsv files from qc_tss_enrichment.sh (default ".")
  --samplesheet  samplesheet.csv; column sample_id is required. condition, replicate, tmm_group and peak_group are
                 copied from the sheet; a column missing from the sheet is filled with blank (metadata only)
  --out          output table (default tss_enrichment.tsv)
Output: tss_enrichment.tsv, one row per sample: sample_id, condition, replicate, tmm_group, peak_group, then every
  column of the per-sample files (n_regions, window_up, window_down, bin_bp, tss_signal, flank_mean, enrichment,
  profile_max, profile_min, max_over_min). The sample comes from the file content, never from the file name.
Exit 3: no input files, malformed input, or a sample_id that is not in the samplesheet.
"""
import argparse
import csv
import os
import sys

META = ["sample_id", "condition", "replicate", "tmm_group", "peak_group"]


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default=".")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", default="tss_enrichment.tsv")
    a = ap.parse_args()

    with open(a.samplesheet, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        ss_rows = list(reader)
    if "sample_id" not in cols:
        die("samplesheet lacks column sample_id")
    # optional metadata columns: blank when the column is not in the sheet
    ss = {r["sample_id"]: r for r in ss_rows}

    files = sorted(f for f in os.listdir(a.input_dir) if f.endswith(".tss_enrichment.tsv"))
    if not files:
        die("no *.tss_enrichment.tsv in %s" % a.input_dir)
    header, rows = None, []
    for f in files:
        with open(os.path.join(a.input_dir, f), newline="", encoding="utf-8") as fh:
            recs = list(csv.DictReader(fh, delimiter="\t"))
        if len(recs) != 1 or "sample_id" not in recs[0]:
            die("%s is not a one-row TSS enrichment file" % f)
        r = recs[0]
        if r["sample_id"] not in ss:
            die("sample %s (from %s) is not in the samplesheet" % (r["sample_id"], f))
        vals = [c for c in r if c != "sample_id"]
        if header is None:
            header = vals
        elif vals != header:
            die("%s has different columns than the other files" % f)
        meta = [(ss[r["sample_id"]].get(c) or "") if c in cols else "" for c in META]
        meta[0] = r["sample_id"]
        rows.append(meta + [r[c] for c in header])
    rows.sort(key=lambda x: x[0])
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(META + header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print("wrote %s (%d samples)" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
