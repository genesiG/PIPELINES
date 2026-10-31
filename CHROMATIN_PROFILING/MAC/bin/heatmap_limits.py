#!/usr/bin/env python3
r"""heatmap_limits.py -- per-sample colour limits of one heatmap, from its computeMatrix matrix.

Called by process HEATMAP_LIMITS (a heatmap whose z_max / y_max entry is 'auto'):
    heatmap_limits.py --matrix NAME.matrix.gz --out-prefix NAME \
        --samples id1,id2 --antibodies ab1,ab2

Arguments:
  --matrix      deepTools computeMatrix .gz of the plotted heatmap (all region groups)
  --out-prefix  NAME (the heatmap id); output NAME_limits.tsv
  --samples     sample ids, comma list, in the order of the matrix sample_labels
  --antibodies  antibody of each sample, comma list, same order as --samples; the antibody string is the
                grouping key: samples with the same string share one limit

Algorithm: per sample and per region group (header group_boundaries), the mean profile over the group's regions
per bin column (NaN ignored), i.e. the lines plotHeatmap draws; profile_max = the maximum over all groups' profiles
(NaN ignored), so no group's profile runs off the plot. Limit per antibody = round_up_2sf(1.05 * the largest
profile_max among that antibody's samples). round_up_2sf rounds up to 2 significant figures (a non-positive or
non-finite input gives 1). A sample whose profile_max is NaN gets limit 1.
The matrix is streamed: a json header line starting with '@' (sample_labels, sample_boundaries, group_boundaries;
a header without group_boundaries is read as one group), then one row per region with 6 BED columns followed by the value columns.

Output (working directory):
  NAME_limits.tsv  header 'sample<TAB>antibody<TAB>profile_max<TAB>limit', one row per sample in matrix order.
                   limit is written with 10 significant digits. The same table is printed to stderr.
Exit 3: matrix missing, empty or not a readable gzip; the header is not a json line starting with '@';
        sample_labels differ from --samples; sample_boundaries do not fit sample_labels; no data rows; a data row
        with a different number of value columns; a value that is not a number; --samples empty or with an empty
        id; the number of --antibodies differs from the number of --samples.
Exit 2: a required argument is missing.
"""
import argparse
import gzip
import json
import math
import os
import sys

import numpy as np


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def split_plain(s):
    if s is None or s.strip() == "":
        return []
    return [t.strip() for t in s.split(",")]


def round_up_2sf(x):
    """Round x up to 2 significant figures; non-positive or non-finite -> 1.0."""
    if x is None or not math.isfinite(x) or x <= 0:
        return 1.0
    e = math.floor(math.log10(x))
    f = 10 ** (1 - e)
    return math.ceil(x * f - 1e-9) / f


def read_profile_max(path, samples):
    """Stream the matrix; return profile_max (float, NaN if none) per sample, in matrix order."""
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            first = fh.readline()
            if not first.startswith("@"):
                die("matrix header line missing (expected a line starting with '@'): %s" % path)
            try:
                header = json.loads(first[1:])
            except ValueError:
                die("matrix header is not valid json: %s" % path)
            labels = header.get("sample_labels")
            bounds = header.get("sample_boundaries")
            if not isinstance(labels, list) or not isinstance(bounds, list) or len(bounds) != len(labels) + 1:
                die("matrix header has no usable sample_labels / sample_boundaries: %s" % path)
            if [str(l) for l in labels] != samples:
                die("sample_labels in the matrix (%s) do not match --samples (%s)"
                    % (", ".join(str(l) for l in labels), ", ".join(samples)))
            width = int(bounds[-1])
            if width <= 0:
                die("matrix has no value columns: %s" % path)

            gb = header.get("group_boundaries")
            if not isinstance(gb, list) or len(gb) < 2:
                gb = None
            ngroups = len(gb) - 1 if gb else 1
            sums = np.zeros((ngroups, width))
            counts = np.zeros((ngroups, width))
            nrows = 0
            g = 0
            for line in fh:
                cols = line.rstrip("\r\n").split("\t")
                nrows += 1
                while gb and g < ngroups - 1 and nrows > int(gb[g + 1]):
                    g += 1
                if len(cols) - 6 != width:
                    die("data row %d has %d value columns; the matrix header needs %d (%s)"
                        % (nrows, len(cols) - 6, width, path))
                try:
                    vals = np.array(cols[6:], dtype=float)
                except ValueError:
                    die("data row %d has a value that is not a number (%s)" % (nrows, path))
                ok = ~np.isnan(vals)
                sums[g] += np.where(ok, vals, 0.0)
                counts[g] += ok
    except (OSError, EOFError) as e:
        die("cannot read matrix %s: %s" % (path, e))

    if nrows == 0:
        die("no data rows in matrix: %s" % path)

    mean = np.full(sums.shape, np.nan)
    has = counts > 0
    mean[has] = sums[has] / counts[has]

    profile_max = []
    for i in range(len(samples)):
        seg = mean[:, int(bounds[i]):int(bounds[i + 1])]
        seg = seg[~np.isnan(seg)]
        profile_max.append(float(seg.max()) if seg.size else float("nan"))
    return profile_max


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--antibodies", required=True)
    a = ap.parse_args()

    if not os.path.isfile(a.matrix) or os.path.getsize(a.matrix) == 0:
        die("matrix missing or empty: %s" % a.matrix)

    samples = split_plain(a.samples)
    antibodies = split_plain(a.antibodies)
    if not samples or any(s == "" for s in samples):
        die("--samples is empty or has an empty id")
    if len(antibodies) != len(samples):
        die("--antibodies has %d entries but --samples has %d" % (len(antibodies), len(samples)))

    profile_max = read_profile_max(a.matrix, samples)

    group_max = {}
    for ab, pm in zip(antibodies, profile_max):
        if math.isfinite(pm):
            group_max[ab] = max(group_max.get(ab, -math.inf), pm)
    group_limit = {ab: round_up_2sf(1.05 * gm) for ab, gm in group_max.items()}
    limits = [("%.10g" % group_limit[ab]) if math.isfinite(pm) else "1"
              for ab, pm in zip(antibodies, profile_max)]

    out_path = a.out_prefix + "_limits.tsv"
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("sample\tantibody\tprofile_max\tlimit\n")
        for s, ab, pm, lim in zip(samples, antibodies, profile_max, limits):
            out.write("%s\t%s\t%.6g\t%s\n" % (s, ab, pm, lim))

    sys.stderr.write("Auto limits for %s:\n" % a.matrix)
    for s, ab, pm, lim in zip(samples, antibodies, profile_max, limits):
        sys.stderr.write("   %-28s %-16s profile_max=%-10.4g limit=%s\n" % (s, ab, pm, lim))


if __name__ == "__main__":
    main()
