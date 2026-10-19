#!/usr/bin/env python3
"""
build_mageck_library.py
========================
Converts a library pool CSV (fixed header: gene,gRNA_ID,sequence) into the
tab-separated <sgRNA_id  sequence  gene> format `mageck count --list-seq`
expects.

Ported from Scripts/step_2_count.py's _build_mageck_library(), made
project-agnostic: column names are fixed (no more GENE_SYMBOL_COL /
ID_COL / TARGET_SEQUENCE_COL index configuration) — every library CSV in
the reusable pipeline must use these exact header names.
"""

import argparse
import csv
import sys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Library pool CSV (gene,gRNA_ID,sequence)")
    ap.add_argument("--output", required=True, help="Output path for the MAGeCK-format TSV")
    args = ap.parse_args()

    with (
        open(args.input, newline="", encoding="utf-8-sig") as fh_in,
        open(args.output, "w", encoding="utf-8") as fh_out,
    ):
        reader = csv.DictReader(fh_in)
        required = {"gene", "gRNA_ID", "sequence"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"ERROR: {args.input} is missing required column(s) {sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        fh_out.write("sgRNA_id\tsequence\tgene\n")
        n = 0
        for i, row in enumerate(reader, start=1):
            gene = (row["gene"] or "").strip().strip("\r")
            gid = (row["gRNA_ID"] or "").strip().strip("\r") or f"{gene}_{i}"
            seq = (row["sequence"] or "").strip().strip("\r")
            if not seq:
                continue
            fh_out.write(f"{gid}\t{seq}\t{gene}\n")
            n += 1

    print(f"Wrote {n} guides to {args.output}")


if __name__ == "__main__":
    main()
