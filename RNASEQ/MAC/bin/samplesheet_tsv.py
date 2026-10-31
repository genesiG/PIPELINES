#!/usr/bin/env python3
"""samplesheet_tsv.py -- the samplesheet as a tab-separated file: kept samples only, the exclude column removed.

Called by SAMPLESHEET_TSV (modules/local/samplesheet_tsv.nf):
    samplesheet_tsv.py --samplesheet samplesheet.csv --keep a,b,c --out samplesheet.tsv

Arguments:
  --samplesheet  the CSV samplesheet (header row required; a UTF-8 byte-order mark is tolerated)
  --keep         comma-separated sample_id values that are processed (the non-excluded rows)
  --out          output TSV
Output: the header row (all columns in sheet order except "exclude"), then the rows whose sample_id is in --keep, in
  the order of the sheet. Tab separated. Metadata columns are carried through unchanged.
Exit 3: missing sheet, no sample_id column, duplicated header names, an empty --keep, or a --keep id that is not in
  the sheet.
"""
import argparse
import csv
import sys


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--keep", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    keep = [k.strip() for k in a.keep.split(",") if k.strip()]
    if not keep:
        die("--keep is empty")
    if len(set(keep)) != len(keep):
        die("duplicate sample ids in --keep")
    keep_set = set(keep)

    try:
        with open(a.samplesheet, encoding="utf-8-sig", newline="") as fh:
            rows = list(csv.reader(fh))
    except OSError as exc:
        die("cannot read samplesheet %s: %s" % (a.samplesheet, exc))
    if not rows:
        die("samplesheet %s is empty" % a.samplesheet)

    header = rows[0]
    names = [h.strip() for h in header]
    if len(set(names)) != len(names):
        die("duplicate header column(s) in %s" % a.samplesheet)
    if "sample_id" not in names:
        die("no sample_id column in %s" % a.samplesheet)
    id_col = names.index("sample_id")
    keep_cols = [i for i in range(len(header)) if names[i] != "exclude"]

    out_rows, seen = [], set()
    for row in rows[1:]:
        if not any(cell.strip() for cell in row):
            continue
        sid = row[id_col].strip() if id_col < len(row) else ""
        if sid in keep_set:
            seen.add(sid)
            out_rows.append([row[i] if i < len(row) else "" for i in keep_cols])

    missing = [k for k in keep if k not in seen]
    if missing:
        die("sample_id(s) in --keep not found in %s: %s" % (a.samplesheet, ", ".join(missing)))

    with open(a.out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([header[i] for i in keep_cols])
        writer.writerows(out_rows)
    print("wrote %s (%d samples)" % (a.out, len(out_rows)))


if __name__ == "__main__":
    main()
