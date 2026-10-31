#!/usr/bin/env python3
"""
library_breakdown.py
=====================
How each sample's reads split across the run's libraries: for every sample,
how many of its reads matched a guide of each library it was counted
against, what share of all its reads that is, and how many matched no
library at all. Most useful when every sample is counted against every
library (no samplesheet `library` column): it shows which pool(s) a sample
actually holds, catches samples that belong to no library (e.g. an empty or
failed sample) and mixtures of pools.

Reads only mageck count outputs, so nothing is re-counted:
  - each library's count table (<prefix>.count.txt: sgRNA, Gene, one column per sample)
  - each library's <prefix>.sample_reads.tsv (bin/mageck_count.sh: sample, fastq,
    reads, mapped, status — status dropped_zero_mapped = the sample had 0 reads
    for that library and was left out of its count table)
  - each library's MAGeCK-format library (sgRNA_id, sequence, gene), so a guide
    sequence present in several libraries is counted once in "Any library".

USAGE:
  python3 library_breakdown.py --libraries=L1,L2 --count-tables="a.count.txt b.count.txt" \
      --sample-reads="a.sample_reads.tsv b.sample_reads.tsv" \
      --mageck-libraries="L1_mageck_library.txt L2_mageck_library.txt" \
      --samples=s1,s2 --project=NAME [--outdir=.]
  --libraries, --count-tables, --sample-reads and --mageck-libraries are in the
  same order (one entry per library). --samples gives the samplesheet order.

OUTPUTS:
  <project>_library_breakdown.csv       long: one row per sample x library (+ rows
                                        "Any library" and "No library guide"):
                                        sample, library, reads, mapped_reads,
                                        pct_of_sample_reads, pct_of_matched_reads, status
  <project>_library_breakdown_wide.csv  one row per sample: reads, % of reads per
                                        library, % any library, % no library guide,
                                        main_library, main_library_pct_of_matched,
                                        n_libraries_ge_5pct_of_matched, shared_guide_reads
"""

import argparse
import csv
import os
import re
import sys

ANY, NONE = "Any library", "No library guide"


def natural_key(name: str):
    """Pool_2 before Pool_10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def split(value: str, sep: str) -> list[str]:
    return [x for x in value.split(sep) if x]


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--libraries", required=True)
    ap.add_argument("--count-tables", required=True)
    ap.add_argument("--sample-reads", required=True)
    ap.add_argument("--mageck-libraries", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    libs = split(args.libraries, ",")
    tables = split(args.count_tables, " ")
    reads_files = split(args.sample_reads, " ")
    lib_files = split(args.mageck_libraries, " ")
    samples = split(args.samples, ",")
    if not (len(libs) == len(tables) == len(reads_files) == len(lib_files)):
        sys.exit(f"ERROR: --libraries ({len(libs)}), --count-tables ({len(tables)}), --sample-reads "
                 f"({len(reads_files)}) and --mageck-libraries ({len(lib_files)}) must have one entry per library.")

    total_reads = {}                       # sample -> input reads (max over libraries; same FASTQ)
    per_lib = {}                           # (sample, lib) -> (mapped, status)
    seq_counts = {}                        # sample -> {sequence: count}
    seq_libs = {}                          # sequence -> set(libraries)
    for lib, table, rfile, lfile in zip(libs, tables, reads_files, lib_files):
        id2seq = {r["sgRNA_id"]: r["sequence"].strip().upper() for r in read_tsv(lfile)}
        for seq in id2seq.values():
            seq_libs.setdefault(seq, set()).add(lib)
        for r in read_tsv(rfile):
            s = r["sample"]
            total_reads[s] = max(total_reads.get(s, 0), int(r["reads"] or 0))
            per_lib[(s, lib)] = (int(r["mapped"] or 0), r["status"])
        for r in read_tsv(table):
            for s, v in r.items():
                if s in ("sgRNA", "Gene"):
                    continue
                seq = id2seq.get(r["sgRNA"])
                if seq is None:
                    continue
                d = seq_counts.setdefault(s, {})
                d[seq] = max(d.get(seq, 0), int(float(v)))   # same read counted in two libraries -> once

    libs = sorted(libs, key=natural_key)   # column/row order only; data is keyed by name
    order = [s for s in samples if s in total_reads] + sorted(s for s in total_reads if s not in samples)
    long_rows, wide_rows = [], []
    for s in order:
        n = total_reads[s]
        counts = seq_counts.get(s, {})
        any_lib = sum(counts.values())
        shared = sum(c for q, c in counts.items() if len(seq_libs.get(q, ())) > 1)
        lib_mapped = {lib: per_lib[(s, lib)] for lib in libs if (s, lib) in per_lib}
        for lib, (m, status) in lib_mapped.items():
            long_rows.append({"sample": s, "library": lib, "reads": n, "mapped_reads": m,
                              "pct_of_sample_reads": round(100 * m / n, 3) if n else 0.0,
                              "pct_of_matched_reads": round(100 * m / any_lib, 3) if any_lib else 0.0,
                              "status": status})
        long_rows.append({"sample": s, "library": ANY, "reads": n, "mapped_reads": any_lib,
                          "pct_of_sample_reads": round(100 * any_lib / n, 3) if n else 0.0,
                          "pct_of_matched_reads": 100.0 if any_lib else 0.0, "status": ""})
        long_rows.append({"sample": s, "library": NONE, "reads": n, "mapped_reads": n - any_lib,
                          "pct_of_sample_reads": round(100 * (n - any_lib) / n, 3) if n else 0.0,
                          "pct_of_matched_reads": "", "status": ""})
        ranked = sorted(lib_mapped.items(), key=lambda kv: kv[1][0], reverse=True)
        main_lib, (main_m, _) = ranked[0] if ranked and ranked[0][1][0] > 0 else ("none", (0, ""))
        row = {"sample": s, "reads": n}
        for lib in libs:
            row[f"pct_{lib}"] = round(100 * lib_mapped[lib][0] / n, 3) if lib in lib_mapped and n else ""
        row.update({"pct_any_library": round(100 * any_lib / n, 3) if n else 0.0,
                    "pct_no_library_guide": round(100 * (n - any_lib) / n, 3) if n else 0.0,
                    "main_library": main_lib,
                    "main_library_pct_of_matched": round(100 * main_m / any_lib, 2) if any_lib else 0.0,
                    "n_libraries_ge_5pct_of_matched": sum(1 for _, (m, _) in lib_mapped.items() if any_lib and m / any_lib >= 0.05),
                    "shared_guide_reads": shared})
        wide_rows.append(row)

    os.makedirs(args.outdir, exist_ok=True)
    for name, rows in ((f"{args.project}_library_breakdown.csv", long_rows),
                       (f"{args.project}_library_breakdown_wide.csv", wide_rows)):
        path = os.path.join(args.outdir, name)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"  Saved: {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
