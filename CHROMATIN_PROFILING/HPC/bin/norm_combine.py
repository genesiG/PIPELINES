#!/usr/bin/env python3
"""norm_combine.py -- combine the per-method scale-factor tables of ONE normalisation run. Standard library only.

CLI contract (called by the NORM_COMBINE module; --key value):

  norm_combine.py --run <name>                       run label (messages only)
                  --primary <method>                 the run's primary method (chipseqspikeinfree|tmm|cpm|spikein)
                  --methods-tsv m=path [m=path ...]  one 7-column table per method that ran (space and/or comma separated)
                  --expected-methods m1,m2,...       methods that were requested for this run (primary + norm_methods +
                                                     spikein when toggled). Each must have a table (--methods-tsv) or a
                                                     reason (--unavailable), otherwise exit 3. A method that is in neither
                                                     list is reported as not_run.
                  [--unavailable m=reason ...]       methods that were not launched, with the reason (space separated; a
                                                     reason may contain spaces if quoted and commas); status unavailable,
                                                     note "UNAVAILABLE: <reason>"
                  --out-prefix <prefix>              output prefix (the run name)
Outputs:
  <prefix>.scale_factors.tsv       BYTE COPY of the primary method's table (schema, path and content as before)
  <prefix>.scale_factors_all.tsv   wide table, tab-separated, one row per library in the primary table's order:
                                   sample_id, primary_method, then for m in chipseqspikeinfree, tmm, cpm, spikein (fixed order):
                                   <m>_raw_factor <m>_bigwig_multiplier <m>_edger_norm_factor <m>_qc_pass <m>_status <m>_note
    <m>_status: ok | qc_fail (finite factors, qc_pass FALSE) | imputed (chipseqspikeinfree only: qc_pass FALSE, raw_factor NA,
                factors finite) | unavailable (no finite factors: note starts "UNAVAILABLE:", or not launched with a reason) |
                not_run (not requested; factor cells NA, qc_pass NA, note "not requested")
Exit codes: 0 ok; 2 usage error; 3 data problem: primary table missing/unavailable or with a missing/NaN/<=0 factor,
  inconsistent sample ids between tables, an expected method with neither table nor reason, unreadable table.
"""

import argparse
import csv
import math
import shutil
import sys

METHODS = ["chipseqspikeinfree", "tmm", "cpm", "spikein"]
COLS = ["sample_id", "method", "raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "note"]
FIELDS = ["raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "status", "note"]


def die(msg, status=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(status)


def num(x):
    try:
        v = float(x)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def read_table(method, path):
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
    except OSError as e:
        die("cannot read the %s table: %s" % (method, e))
    if not rows or rows[0] != COLS:
        die("%s: header is %s, expected %s" % (path, rows[0] if rows else None, COLS))
    body = rows[1:]
    if not body or any(len(r) != 7 for r in body):
        die("%s: no rows or a row without 7 columns" % path)
    ids = [r[0] for r in body]
    if len(set(ids)) != len(ids):
        die("%s: duplicated sample_id" % path)
    for r in body:
        if r[1] != method:
            die("%s: method column is '%s', expected '%s'" % (path, r[1], method))
    return {r[0]: r for r in body}, ids


def status_of(method, r):
    bw, nf = num(r[3]), num(r[4])
    if r[6].startswith("UNAVAILABLE:") or bw is None or nf is None:
        return "unavailable"
    if r[5] == "TRUE":
        return "ok"
    if method == "chipseqspikeinfree" and num(r[2]) is None:
        return "imputed"
    return "qc_fail"


def split_items(values):
    out = []
    for v in values:
        out.extend(x.strip() for x in v.split(",") if x.strip())
    return out


def main():
    ap = argparse.ArgumentParser(description="Combine per-method scale-factor tables of one run.")
    ap.add_argument("--run", required=True)
    ap.add_argument("--primary", required=True)
    ap.add_argument("--methods-tsv", nargs="+", required=True)
    ap.add_argument("--expected-methods", nargs="+", required=True)
    ap.add_argument("--unavailable", nargs="*", default=[])
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    if a.primary not in METHODS:
        die("unknown primary method '%s' (allowed: %s)" % (a.primary, ", ".join(METHODS)), 2)
    paths = {}
    for item in split_items(a.methods_tsv):
        if "=" not in item:
            die("--methods-tsv items must be method=path, got '%s'" % item, 2)
        m, p = item.split("=", 1)
        if m not in METHODS:
            die("unknown method '%s' in --methods-tsv" % m, 2)
        paths[m] = p
    expected = split_items(a.expected_methods)
    for m in expected:
        if m not in METHODS:
            die("unknown method '%s' in --expected-methods" % m, 2)
    reasons = {}
    for item in a.unavailable:   # not split on commas: a reason is free text
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            die("--unavailable items must be method=reason, got '%s'" % item, 2)
        m, why = item.split("=", 1)
        if m not in METHODS:
            die("unknown method '%s' in --unavailable" % m, 2)
        reasons[m] = why.strip() or "not available"

    # ---- primary check
    if a.primary not in paths:
        if a.primary in reasons:
            die("run %s: the primary method %s is unavailable: %s" % (a.run, a.primary, reasons[a.primary]))
        die("run %s: no scale-factor table for the primary method %s" % (a.run, a.primary))
    tables = {}
    for m, p in paths.items():
        tables[m], ids = read_table(m, p)
        if m == a.primary:
            order = ids   # run order = the primary table's row order
    for m, t in tables.items():
        if set(t) != set(order):
            die("run %s: the %s table has different sample ids than the primary table" % (a.run, m))
    prim = tables[a.primary]
    bad = []
    for sid in order:
        r = prim[sid]
        if r[6].startswith("UNAVAILABLE:"):
            die("run %s: the primary method %s is unavailable (%s)" % (a.run, a.primary, r[6]))
        for col, name in ((3, "bigwig_multiplier"), (4, "edger_norm_factor")):
            v = num(r[col])
            if v is None or v <= 0:
                bad.append("%s %s=%s" % (sid, name, r[col]))
    if bad:
        die("run %s: the primary method %s has missing, NaN or non-positive scale factors for: %s. "
            "Fix the input or choose another primary method." % (a.run, a.primary, "; ".join(bad)))
    for m in expected:
        if m not in paths and m not in reasons:
            die("run %s: method %s was expected but has neither a table nor an unavailable reason" % (a.run, m))

    # ---- outputs
    shutil.copyfile(paths[a.primary], a.out_prefix + ".scale_factors.tsv")
    header = ["sample_id", "primary_method"] + ["%s_%s" % (m, f) for m in METHODS for f in FIELDS]
    with open(a.out_prefix + ".scale_factors_all.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(header) + "\n")
        for sid in order:
            cells = [sid, a.primary]
            for m in METHODS:
                if m in tables:
                    r = tables[m][sid]
                    cells += [r[2], r[3], r[4], r[5], status_of(m, r), r[6]]
                elif m in reasons:
                    cells += ["NA", "NA", "NA", "FALSE", "unavailable", "UNAVAILABLE: " + reasons[m]]
                else:
                    cells += ["NA", "NA", "NA", "NA", "not_run", "not requested"]
            fh.write("\t".join(cells) + "\n")
    print("Wrote %s.scale_factors.tsv (copy of %s) and %s.scale_factors_all.tsv (%d libraries; methods run: %s)"
          % (a.out_prefix, a.primary, a.out_prefix, len(order), ", ".join(m for m in METHODS if m in tables) or "none"))


if __name__ == "__main__":
    main()
