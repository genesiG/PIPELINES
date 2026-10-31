#!/usr/bin/env python3
"""heatmap_cluster.py -- k-means clustering of the regions of a heatmap matrix.

Called by HEATMAP_CLUSTER (modules/local/heatmap_cluster.nf, container cutandrun-py_1.1.sif):
    heatmap_cluster.py --matrix NAME.matrix.gz --regions regions.bed --samples id1,id2 \
        --k 2-6 --seed 42 --out-prefix NAME

Arguments (--key value):
  --matrix        deepTools computeMatrix output (.gz or plain text; first character '@' followed
                  by the JSON header), or a plain numeric TSV without header (gz or plain).
                  The tsv format is detected by the first character not being '@'.
  --regions       original region BED. Its lines are copied into the cluster BEDs. Optional for
                  deepTools matrices; when given, each matrix row is matched to the BED line with the
                  same chrom, start and end (first three fields), not by position, because deepTools
                  may drop or reorder regions (see Matrix format). REQUIRED for the tsv format, which
                  needs the same row order and row count as the matrix.
  --samples       deepTools matrices only: comma-separated sample_labels whose columns are
                  clustered. Empty = all samples. Ignored with a warning for the tsv format.
  --k             "3" (fixed k) or "2-6" (range; the gap statistic chooses k). Minimum k is 2.
                  A range with equal ends (e.g. "4-4") is treated as a fixed k.
  --seed          non-negative integer. Seeds KMeans, the silhouette subsample and the gap
                  reference datasets.
  --out-prefix    NAME = the heatmap id. Prefix of all outputs; written to the 'heatmap' column.
  --n-init        KMeans restarts for the observed data (default 100).
  --max-iter      KMeans max iterations (default 500).
  --gap-refs      uniform reference datasets for the gap statistic (default 20, must be >= 2;
                  used only for a range of k).
  --sil-sample    maximum number of regions in the silhouette subsample (default 5000).
  --threads       parallel workers for the gap reference runs (joblib; 1 = serial, default 1).

Matrix format (deepTools)
  Line 1: '@' + JSON with keys sample_labels, sample_boundaries (n_samples+1 offsets into the
  value columns, starting at 0), group_labels and group_boundaries. Data rows: 6 metadata columns
  (chrom, start, end, name, score, strand), then the values ('nan' allowed). Only the columns of
  the requested samples are read; the metadata is kept for the BED6 output.
  With --regions, row i of the matrix is matched to the next unused --regions line with the same
  (chrom, start, end) key (in BED order). Matrix rows fewer than the BED lines give a WARNING (regions
  dropped by deepTools); a matrix row without a match is an error (exit 3).

Outputs (NAME = --out-prefix)
  NAME.cluster_assignments.csv  one integer per matrix row (1-indexed, no header, matrix row
                                order) for the chosen k, after relabelling.
  NAME.cluster_k.tsv            one row per k evaluated. Columns: heatmap k silhouette gap gap_se
                                chosen. gap and gap_se are NA for a fixed k; silhouette is NA when
                                the subsample holds fewer than 2 labels. chosen is TRUE or FALSE.
  NAME.clusters.tsv             chosen k only. Columns: heatmap k cluster n_regions mean_signal.
  NAME.cluster<c>.bed           c = 1..chosen k. With --regions: the original lines in matrix order
                                (comment, track and browser lines and blank lines are dropped).
                                Without --regions: BED6 from the deepTools row metadata.
  Numbers are written with 6 decimals.

Exit codes
  0  ok
  2  usage error (missing or malformed arguments, k spec, tsv input without --regions)
  3  data problem: empty or all-NaN matrix, k > number of regions, tsv format with --regions line
     count not equal to the matrix row count, deepTools matrix row with no matching --regions line
     (key chrom, start, end), sample name not in sample_labels, empty cluster in the chosen k,
     unreadable or malformed input.

Method
  Matrix: the values of the requested columns are read as float32; NaN is set to 0. This matrix is
  the clustering matrix (also used for the cluster mean_signal and the relabelling).
  KMeans: sklearn KMeans(n_clusters=k, init='k-means++', n_init=--n-init, max_iter=--max-iter,
  random_state=--seed, algorithm='lloyd') for every k in the spec.
  Silhouette: Euclidean silhouette on a seeded subsample of at most --sil-sample regions
  (np.random.default_rng(seed)).
  Gap statistic (Tibshirani et al. 2001): only for a range of k. Reference datasets are uniform
  within the bounding box of the clustering matrix. Each reference is clustered with
  KMeans(n_init=max(1, n_init // 5)) and random_state = seed + i*97 + k*13 (i = 0..B-1, B = gap
  refs). gap(k) = mean(log W_ref) - log W_obs; se(k) = sd(log W_ref, ddof=1) * sqrt(1 + 1/B).
  Chosen k (Tibshirani 1-SE rule): the smallest k with gap(k) >= gap(k+1) - se(k+1). If no k
  satisfies it, the k with the largest gap.
  Relabelling: the clusters of the chosen k are renumbered 1..k in decreasing order of the mean of
  the clustering matrix over their rows. Ties: the larger cluster first, then the lower original
  label. The same labels are written to the assignments, the clusters table and the BEDs.

Why deepTools --kmeans is not used
  deepTools 3.5.6 (heatmapper.py:1279) calls scipy.cluster.vq.kmeans without a seed, so its
  cluster assignments change from run to run. It also reports no silhouette or gap statistic.
  This script seeds every random step, so the same input gives the same clusters.
  Port of the clustering logic of the CBX2 step_8e_clustering.py (run_clustering) and of the
  matrix extraction idea of CBX2 utils/extract_cluster_matrix.R.
"""
import argparse
import gzip
import itertools
import json
import os
import re
import sys

import numpy as np
from joblib import Parallel, delayed
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

DT_KEYS = ["sample_labels", "sample_boundaries", "group_labels", "group_boundaries"]
META_COLS = 6
BED_SKIP = ("#", "track", "browser")


def die(msg, code=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(code)


def open_text(path):
    """Open a gzip or plain text file (detected by the gzip magic bytes 1f 8b)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def peek_first_line(path):
    if not os.path.isfile(path):
        die("matrix not found: %s" % path)
    with open_text(path) as fh:
        first = fh.readline()
    if not first:
        die("matrix is empty: %s" % path)
    return first


def parse_k(spec):
    """'3' -> ([3], False); '2-6' -> ([2..6], True); None if malformed or kmin < 2 or kmin > kmax."""
    m = re.fullmatch(r"\s*(\d+)(?:-(\d+))?\s*", spec)
    if not m:
        return None
    kmin = int(m.group(1))
    kmax = int(m.group(2)) if m.group(2) is not None else kmin
    if not (2 <= kmin <= kmax):
        return None
    return list(range(kmin, kmax + 1)), kmin != kmax


def row_to_float(parts, cols, where):
    try:
        return np.array([float(parts[c]) for c in cols], dtype=np.float32)
    except ValueError:
        die("non-numeric value at %s" % where)


def read_deeptools(fh, first, samples):
    try:
        hdr = json.loads(first[1:])
    except json.JSONDecodeError as exc:
        die("cannot parse the deepTools header line: %s" % exc)
    for key in DT_KEYS:
        if key not in hdr:
            die("deepTools header lacks the key '%s'" % key)
    labels = hdr["sample_labels"]
    bounds = hdr["sample_boundaries"]
    if len(bounds) != len(labels) + 1 or bounds[0] != 0 or any(
            bounds[i + 1] < bounds[i] for i in range(len(labels))):
        die("sample_boundaries are inconsistent with sample_labels")
    n_val = bounds[-1]

    if samples:
        sel = []
        for name in samples:
            hits = [i for i, lab in enumerate(labels) if lab == name]
            if not hits:
                die("sample '%s' not in the matrix sample_labels (%s)" % (name, ", ".join(labels)))
            if len(hits) > 1 or hits[0] in sel:
                die("sample '%s' is duplicated" % name)
            sel.append(hits[0])
    else:
        sel = list(range(len(labels)))

    value_cols = []
    for j in sel:
        value_cols.extend(range(bounds[j], bounds[j + 1]))
    cols = [META_COLS + c for c in value_cols]
    if not cols:
        die("no value columns selected")
    sel_labels = [labels[j] for j in sel]

    n_fields = META_COLS + n_val
    rows, meta = [], []
    for lineno, line in enumerate(fh, start=2):
        line = line.rstrip("\r\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != n_fields:
            die("matrix line %d has %d fields, expected %d" % (lineno, len(parts), n_fields))
        meta.append(parts[:META_COLS])
        rows.append(row_to_float(parts, cols, "matrix line %d" % lineno))
    return rows, meta, sel_labels


def read_tsv(fh, first):
    rows, ncol = [], None
    for lineno, line in enumerate(itertools.chain([first], fh), start=1):
        line = line.rstrip("\r\n")
        if not line:
            continue
        parts = line.split("\t")
        if ncol is None:
            ncol = len(parts)
        elif len(parts) != ncol:
            die("matrix line %d has %d columns, expected %d" % (lineno, len(parts), ncol))
        rows.append(row_to_float(parts, range(ncol), "matrix line %d" % lineno))
    return rows, None, None


def load_matrix(path, is_dt, samples):
    """Return (X float32 with NaN -> 0, row metadata or None, sample labels or None)."""
    with open_text(path) as fh:
        first = fh.readline()
        if is_dt:
            rows, meta, labels = read_deeptools(fh, first, samples)
        else:
            rows, meta, labels = read_tsv(fh, first)
    if not rows:
        die("matrix has no data rows: %s" % path)
    X = np.vstack(rows)
    del rows
    nan = np.isnan(X)
    if nan.all():
        die("matrix contains only NaN values: %s" % path)
    X[nan] = 0
    if not np.isfinite(X).all():
        die("matrix contains infinite values: %s" % path)
    return X, meta, labels


def ref_log_inertia(shape, col_min, col_max, k, seed, n_init, max_iter):
    """log(inertia) of KMeans on one uniform reference dataset in the bounding box."""
    rng = np.random.default_rng(seed)
    X_ref = rng.uniform(low=col_min, high=col_max, size=shape).astype(np.float32)
    km = KMeans(
        n_clusters=k,
        init="k-means++",
        n_init=n_init,
        max_iter=max_iter,
        random_state=seed,
        algorithm="lloyd",
    )
    km.fit(X_ref)
    return float(np.log(km.inertia_))


def relabel(X, labels, k):
    """Renumber clusters 1..k by decreasing mean of X over the rows; ties: larger, then lower label."""
    sizes = np.bincount(labels, minlength=k)
    means = np.array([
        X[labels == c].mean(dtype=np.float64) if sizes[c] > 0 else -np.inf
        for c in range(k)
    ])
    order = sorted(range(k), key=lambda c: (-means[c], -sizes[c], c))
    lut = np.zeros(k, dtype=np.int64)
    for new_label, c in enumerate(order, start=1):
        lut[c] = new_label
    return lut[labels]


def fnum(v):
    if v is None or not np.isfinite(v):
        return "NA"
    return "%.6f" % v


def bed_key(fields, where):
    """(chrom, int start, int end) from the first three fields of a BED line or a matrix row."""
    try:
        return (fields[0], int(fields[1]), int(fields[2]))
    except (IndexError, ValueError):
        die("cannot read chrom, start and end at %s" % where)


def match_bed_lines(bed_lines, meta):
    """Return the BED lines in matrix row order, matched by key (chrom, start, end).

    Each matrix row takes the next unused BED line with the same key (in BED order). BED lines
    without a matrix row (regions dropped by deepTools) are skipped with a WARNING; a matrix row
    without an unused BED line is an error.
    """
    index = {}
    for i, s in enumerate(bed_lines):
        index.setdefault(bed_key(s.split("\t"), "region %d of --regions" % (i + 1)), []).append(i)
    used = {}
    matched = []
    for r, m in enumerate(meta):
        key = bed_key(m, "matrix row %d" % (r + 1))
        cand = index.get(key, [])
        j = used.get(key, 0)
        if j >= len(cand):
            die("matrix row %d (%s:%d-%d) has no unused --regions line with the same chrom, start and end"
                % (r + 1, key[0], key[1], key[2]))
        matched.append(cand[j])
        used[key] = j + 1
    n_dropped = len(bed_lines) - len(matched)
    if n_dropped > 0:
        sys.stderr.write("WARNING: %d of %d --regions lines have no matrix row (dropped by deepTools) "
                         "and are not in the cluster BEDs\n" % (n_dropped, len(bed_lines)))
    return [bed_lines[i] for i in matched]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--regions", default=None)
    ap.add_argument("--samples", default="")
    ap.add_argument("--k", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-init", type=int, default=100)
    ap.add_argument("--max-iter", type=int, default=500)
    ap.add_argument("--gap-refs", type=int, default=20)
    ap.add_argument("--sil-sample", type=int, default=5000)
    ap.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()

    parsed = parse_k(a.k)
    if parsed is None:
        ap.error("--k must be 'K' or 'KMIN-KMAX' with 2 <= KMIN <= KMAX, got '%s'" % a.k)
    k_values, is_range = parsed
    if a.seed < 0:
        ap.error("--seed must be >= 0")
    if a.n_init < 1 or a.max_iter < 1:
        ap.error("--n-init and --max-iter must be >= 1")
    if a.sil_sample < 2:
        ap.error("--sil-sample must be >= 2")
    if a.threads < 1:
        ap.error("--threads must be >= 1")
    if is_range and a.gap_refs < 2:
        ap.error("--gap-refs must be >= 2 for a range of k")

    first = peek_first_line(a.matrix)
    is_dt = first.startswith("@")
    if not is_dt:
        if not a.regions:
            ap.error("--regions is required for a tsv-format matrix (no deepTools header)")
        if a.samples.strip():
            sys.stderr.write("WARNING: --samples is ignored for a tsv-format matrix\n")
    samples = [s.strip() for s in a.samples.split(",") if s.strip()] if is_dt else []

    X, meta, labels = load_matrix(a.matrix, is_dt, samples)
    n, p = X.shape
    fmt_name = "deepTools" if is_dt else "tsv"
    print("=== heatmap_cluster.py: %s ===" % a.out_prefix)
    print("matrix: %s (%s format), %d regions x %d columns" % (a.matrix, fmt_name, n, p))
    if is_dt:
        print("clustering samples: %s" % ", ".join(labels))

    if max(k_values) > n:
        die("k=%d exceeds the number of regions (%d)" % (max(k_values), n))

    bed_lines = None
    if a.regions:
        if not os.path.isfile(a.regions):
            die("regions BED not found: %s" % a.regions)
        bed_lines = []
        with open(a.regions, encoding="utf-8") as fh:
            for line in fh:
                s = line.rstrip("\r\n")
                if not s.strip() or s.startswith(BED_SKIP):
                    continue
                bed_lines.append(s)
        if not is_dt:
            if len(bed_lines) != n:
                die("--regions has %d region lines but the matrix has %d rows" % (len(bed_lines), n))
        else:
            bed_lines = match_bed_lines(bed_lines, meta)

    # Silhouette subsample (drawn once, as in the legacy script)
    rng = np.random.default_rng(a.seed)
    if n > a.sil_sample:
        sil_idx = rng.choice(n, size=a.sil_sample, replace=False)
        print("silhouette subsample: %d / %d regions" % (a.sil_sample, n))
    else:
        sil_idx = np.arange(n)

    print("KMeans: init=k-means++ n_init=%d max_iter=%d seed=%d" % (a.n_init, a.max_iter, a.seed))
    print("\n--- Fitting KMeans ---")
    inertia = {}
    sil = {}
    labels_by_k = {}
    for k in k_values:
        km = KMeans(
            n_clusters=k,
            init="k-means++",
            n_init=a.n_init,
            max_iter=a.max_iter,
            random_state=a.seed,
            algorithm="lloyd",
        )
        km.fit(X)
        lab = km.labels_
        labels_by_k[k] = lab
        inertia[k] = float(km.inertia_)

        lb_sub = lab[sil_idx]
        n_lab = len(np.unique(lb_sub))
        if n_lab < 2 or n_lab >= len(sil_idx):
            sil[k] = float("nan")
        else:
            sil[k] = float(silhouette_score(X[sil_idx], lb_sub, metric="euclidean"))

        sizes = np.bincount(lab, minlength=k)
        size_txt = " | ".join("C%d: %d" % (c + 1, sizes[c]) for c in range(k))
        print("  k=%d  inertia=%.1f  silhouette=%s  [%s]" % (k, inertia[k], fnum(sil[k]), size_txt))

    gap = {k: None for k in k_values}
    gap_se = {k: None for k in k_values}
    if is_range:
        print("\n--- Gap statistic (%d uniform references, threads=%d) ---" % (a.gap_refs, a.threads))
        col_min = X.min(axis=0)
        col_max = X.max(axis=0)
        ref_n_init = max(1, a.n_init // 5)
        for k in k_values:
            seeds = [a.seed + i * 97 + k * 13 for i in range(a.gap_refs)]
            args_k = [(X.shape, col_min, col_max, k, s, ref_n_init, a.max_iter) for s in seeds]
            if a.threads == 1:
                log_refs = [ref_log_inertia(*t) for t in args_k]
            else:
                log_refs = Parallel(n_jobs=a.threads)(delayed(ref_log_inertia)(*t) for t in args_k)
            log_refs = np.asarray(log_refs, dtype=np.float64)
            gap[k] = float(log_refs.mean() - np.log(inertia[k]))
            gap_se[k] = float(log_refs.std(ddof=1) * np.sqrt(1 + 1 / a.gap_refs))
            print("  k=%d  gap=%.4f  SE=%.4f" % (k, gap[k], gap_se[k]))

        # Tibshirani 1-SE rule: smallest k with gap(k) >= gap(k+1) - se(k+1); else argmax gap
        chosen = k_values[0]
        for i in range(len(k_values) - 1):
            if gap[k_values[i]] >= gap[k_values[i + 1]] - gap_se[k_values[i + 1]]:
                chosen = k_values[i]
                break
        else:
            chosen = k_values[int(np.argmax([gap[k] for k in k_values]))]
        print("\nChosen k = %d (gap statistic, 1-SE rule)" % chosen)
    else:
        chosen = k_values[0]
        print("\nChosen k = %d (fixed)" % chosen)

    new_labels = relabel(X, labels_by_k[chosen], chosen)
    cl_sizes = np.bincount(new_labels - 1, minlength=chosen)
    if (cl_sizes == 0).any():
        die("chosen k=%d has an empty cluster (sizes %s)" % (chosen, cl_sizes.tolist()))

    prefix = a.out_prefix

    # --- NAME.cluster_assignments.csv ---
    with open(prefix + ".cluster_assignments.csv", "w", encoding="utf-8") as fh:
        fh.write("".join("%d\n" % v for v in new_labels))

    # --- NAME.cluster_k.tsv ---
    with open(prefix + ".cluster_k.tsv", "w", encoding="utf-8") as fh:
        fh.write("heatmap\tk\tsilhouette\tgap\tgap_se\tchosen\n")
        for k in k_values:
            fh.write("%s\t%d\t%s\t%s\t%s\t%s\n" % (
                prefix, k, fnum(sil[k]), fnum(gap[k]), fnum(gap_se[k]),
                "TRUE" if k == chosen else "FALSE"))

    # --- NAME.clusters.tsv ---
    with open(prefix + ".clusters.tsv", "w", encoding="utf-8") as fh:
        fh.write("heatmap\tk\tcluster\tn_regions\tmean_signal\n")
        for c in range(1, chosen + 1):
            sel = new_labels == c
            fh.write("%s\t%d\t%d\t%d\t%s\n" % (
                prefix, chosen, c, int(sel.sum()), fnum(X[sel].mean(dtype=np.float64))))

    # --- NAME.cluster<c>.bed ---
    for c in range(1, chosen + 1):
        idx = np.flatnonzero(new_labels == c)
        if bed_lines is not None:
            out_lines = [bed_lines[int(i)] for i in idx]
        else:
            out_lines = ["\t".join(meta[int(i)]) for i in idx]
        with open("%s.cluster%d.bed" % (prefix, c), "w", encoding="utf-8") as fh:
            fh.write("\n".join(out_lines) + "\n")

    print("\nCluster sizes (k=%d, relabelled): %s" % (chosen, " | ".join(
        "C%d: %d" % (c + 1, cl_sizes[c]) for c in range(chosen))))
    print("wrote %s.cluster_assignments.csv, %s.cluster_k.tsv, %s.clusters.tsv, %s.cluster1.bed ... %s.cluster%d.bed"
          % (prefix, prefix, prefix, prefix, prefix, chosen))


if __name__ == "__main__":
    main()
