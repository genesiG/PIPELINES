#!/usr/bin/env python3
"""qc_library_complexity_table.py -- merge the per-sample library complexity results into library_complexity.tsv.

Called by QC_LIBCOMPLEXITY_TABLE (modules/local/qc_libcomplexity_table.nf):
    qc_library_complexity_table.py --input-dir . --samplesheet samplesheet.csv --out library_complexity.tsv

Arguments:
  --input-dir    directory with the per-sample files (default "."):
                   <S>.complexity.tsv       NRF/PBC1/PBC2 from fragment coordinates (LIBCOMPLEXITY_FRAGS)
                   <S>.elc.txt              Picard EstimateLibraryComplexity metrics (PICARD_ELC; optional)
                   <S>.markdup_metrics.txt  Picard MarkDuplicates metrics (optional)
  --samplesheet  samplesheet.csv (columns sample_id, antibody, condition, replicate, group); metadata only from here
  --out          output table (default library_complexity.tsv)
Output: library_complexity.tsv, one row per sample that has an <S>.complexity.tsv:
  sample_id, antibody, condition, replicate, group, total_fragments, distinct_fragments, once, twice, NRF, PBC1,
  PBC2, elc_read_pairs_examined, elc_pct_duplication, elc_estimated_library_size, markdup_read_pairs_examined,
  markdup_pct_duplication, markdup_estimated_library_size  (NA where a file or value is absent)
Exit 3: no <S>.complexity.tsv, or one whose sample_id is not in the samplesheet.
"""
import argparse
import csv
import os
import sys

META = ["sample_id", "antibody", "condition", "replicate", "group"]
COMPLEXITY = ["total_fragments", "distinct_fragments", "once", "twice", "NRF", "PBC1", "PBC2"]


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def read_samplesheet(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or any(c not in rows[0] for c in META):
        die("samplesheet %s lacks one of the columns %s" % (path, ", ".join(META)))
    return {r["sample_id"]: r for r in rows}


def picard_metrics(path):
    """First metrics row of a Picard metrics file as a dict, or {} (absent / header only / unparsable)."""
    if not os.path.exists(path):
        return {}
    lines = [l.rstrip("\n") for l in open(path, encoding="utf-8", errors="replace")]
    for i, l in enumerate(lines):
        if "READ_PAIRS_EXAMINED" in l and i + 1 < len(lines) and lines[i + 1].strip():
            hdr, val = l.split("\t"), lines[i + 1].split("\t")
            if len(val) >= len(hdr):
                return dict(zip(hdr, val))
    return {}


def num(d, key, scale=1.0):
    try:
        return "%.6g" % (float(d[key]) * scale)
    except (KeyError, ValueError):
        return "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", default=".")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", default="library_complexity.tsv")
    a = ap.parse_args()

    ss = read_samplesheet(a.samplesheet)
    files = sorted(f for f in os.listdir(a.input_dir) if f.endswith(".complexity.tsv"))
    if not files:
        die("no *.complexity.tsv in %s" % a.input_dir)

    out = []
    for f in files:
        with open(os.path.join(a.input_dir, f), encoding="utf-8") as fh:
            recs = list(csv.DictReader(fh, delimiter="\t"))
        if len(recs) != 1 or "sample_id" not in recs[0]:
            die("%s must hold a header and exactly one row" % f)
        r = recs[0]
        sid = r["sample_id"]
        if sid not in ss:
            die("sample %s (from %s) is not in the samplesheet" % (sid, f))
        elc = picard_metrics(os.path.join(a.input_dir, sid + ".elc.txt"))
        md = picard_metrics(os.path.join(a.input_dir, sid + ".markdup_metrics.txt"))
        row = [ss[sid][c] for c in META] + [r.get(c, "NA") for c in COMPLEXITY]
        row += [num(elc, "READ_PAIRS_EXAMINED"), num(elc, "PERCENT_DUPLICATION", 100), num(elc, "ESTIMATED_LIBRARY_SIZE"),
                num(md, "READ_PAIRS_EXAMINED"), num(md, "PERCENT_DUPLICATION", 100), num(md, "ESTIMATED_LIBRARY_SIZE")]
        out.append(row)

    header = META + COMPLEXITY + ["elc_read_pairs_examined", "elc_pct_duplication", "elc_estimated_library_size",
                                  "markdup_read_pairs_examined", "markdup_pct_duplication", "markdup_estimated_library_size"]
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in out:
            fh.write("\t".join(row) + "\n")
    print("wrote %s (%d samples)" % (a.out, len(out)))


if __name__ == "__main__":
    main()
