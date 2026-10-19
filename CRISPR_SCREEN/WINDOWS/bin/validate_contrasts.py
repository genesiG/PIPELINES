#!/usr/bin/env python3
"""
validate_contrasts.py
=======================
Resolves each (treatment, control) row in contrasts.csv to the MAGeCK count
table it belongs to and to the actual sample columns in that table, and
fails fast — before any mageck test job is submitted — if it can't.

A run can have several count tables (one per library, when samplesheet.csv
assigns samples to libraries). Counts from different libraries can't be
compared, so each contrast must resolve, both sides, within exactly one
table:
  - a side that matches no column in any table            -> error
  - treatment and control only found in different tables  -> error
  - both sides found together in more than one table      -> error (ambiguous)

Matching is the legacy Scripts/step_3_mageck_test.py rule: regex search,
case-insensitive, trailing "_column\\d+" stripped, --exclude patterns removed.
Because it's a search, "Pool_1" also matches "Pool_10"; anchor a pattern
("^Pool_1$") to match one sample exactly.

Errors exit with status 3 (a data problem: the pipeline doesn't retry it).

Output CSV columns: count_table, name, treatment, control, treat_cols, control_cols
(count_table = file name of the table the contrast was resolved in; name =
<treatment>_vs_<control> with regex characters removed, used for the
contrast's output directory and file names, so it must be unique).
"""

import argparse
import csv
import os
import re
import sys


def fail(message: str) -> None:
    print(message, file=sys.stderr)
    sys.exit(3)


def resolve(pattern: str, header_cols: list[str], exclude: list[str]) -> list[str]:
    cols = [
        c for c in header_cols
        if re.search(pattern, c, re.IGNORECASE)
        and not any(re.search(ex, c, re.IGNORECASE) for ex in exclude)
    ]
    return [re.sub(r"_column\d+", "", c) for c in cols]


def safe_name(pattern: str) -> str:
    """'^Pool_1$' -> 'Pool_1': a pattern made fit for file and directory names."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", pattern.replace("^", "").replace("$", "")).strip("_")


def sample_columns(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as fh:
        header_cols = fh.readline().rstrip("\n").split("\t")
    return [c for c in header_cols if c not in ("sgRNA", "Gene")]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count-table", required=True, nargs="+", help="one or more MAGeCK count tables")
    ap.add_argument("--contrasts", required=True, help="CSV: treatment,control")
    ap.add_argument("--output", required=True)
    ap.add_argument("--exclude", default="", help="comma-separated regex patterns to exclude from matches")
    args = ap.parse_args()

    exclude = [p for p in args.exclude.split(",") if p.strip()]
    tables = {os.path.basename(p): sample_columns(p) for p in sorted(args.count_table)}
    available = "\n".join(f"  {name}: {', '.join(cols)}" for name, cols in tables.items())

    with open(args.contrasts, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        missing = {"treatment", "control"} - set(reader.fieldnames or [])
        if missing:
            fail(f"ERROR: {args.contrasts} is missing required column(s) {sorted(missing)}.")
        rows = [r for r in reader if (r.get("treatment") or "").strip() and (r.get("control") or "").strip()]

    if not rows:
        fail(f"ERROR: {args.contrasts} has no (treatment, control) rows to test.\n"
             f"Sample columns available, per count table:\n{available}")

    errors = []
    resolved = []
    for row in rows:
        treat, ctrl = row["treatment"].strip(), row["control"].strip()
        if treat.lower() == ctrl.lower():
            continue
        per_table = {
            name: (resolve(treat, cols, exclude), resolve(ctrl, cols, exclude))
            for name, cols in tables.items()
        }
        both = [name for name, (t, c) in per_table.items() if t and c]
        if len(both) == 1:
            name = both[0]
            t_cols, c_cols = per_table[name]
            resolved.append({
                "count_table": name,
                "name": f"{safe_name(treat)}_vs_{safe_name(ctrl)}",
                "treatment": treat,
                "control": ctrl,
                "treat_cols": ",".join(t_cols),
                "control_cols": ",".join(c_cols),
            })
            continue
        if len(both) > 1:
            errors.append(f"  '{treat}' vs '{ctrl}': both sides match in more than one count table "
                          f"({', '.join(both)}) — make the patterns specific to one library")
            continue
        t_in = [name for name, (t, _) in per_table.items() if t]
        c_in = [name for name, (_, c) in per_table.items() if c]
        if not t_in:
            errors.append(f"  treatment '{treat}' matched no sample column in any count table")
        if not c_in:
            errors.append(f"  control   '{ctrl}' matched no sample column in any count table")
        if t_in and c_in:
            errors.append(f"  '{treat}' vs '{ctrl}': treatment is in {', '.join(t_in)} but control is in "
                          f"{', '.join(c_in)} — samples screened with different libraries can't be compared")

    names = [r["name"] for r in resolved]
    for dup in sorted({n for n in names if names.count(n) > 1}):
        errors.append(f"  several contrasts would be named '{dup}' (their output folder) — make their patterns differ")

    if errors:
        fail(f"ERROR: {len(errors)} problem(s) in {args.contrasts}:\n"
             + "\n".join(errors)
             + f"\n\nSample columns available, per count table:\n{available}\n"
             "Fix contrasts.csv (or the samplesheet's `library` column) and re-run.")

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["count_table", "name", "treatment", "control", "treat_cols", "control_cols"])
        writer.writeheader()
        writer.writerows(resolved)

    print(f"Validated {len(resolved)} contrast(s) across {len(tables)} count table(s). Wrote {args.output}")


if __name__ == "__main__":
    main()
