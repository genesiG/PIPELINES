#!/usr/bin/env python3
"""qc_tss_enrichment.py -- TSS enrichment of one normalised bigWig (numpy + pyBigWig, no deepTools).

Called by TSS_ENRICHMENT (modules/local/tss_enrichment.nf):
    qc_tss_enrichment.py --sample S --bigwig S.<method>.<units>.bw --gtf genes.gtf --threads N
    qc_tss_enrichment.py --sample S --bigwig S.<method>.<units>.bw --tss-bed tss.bed --threads N   (fallback, no GTF)

Arguments:
  --sample        sample id (required)
  --bigwig        bigWig of the sample (required)
  --gtf           GTF (plain or .gz); TSS = transcript start ('+') or transcript end ('-') of each 'transcript'
                  feature, strand from GTF column 7; for a transcript_id without a 'transcript' feature, min exon
                  start / max exon end. TSSs are deduplicated by (chrom, TSS, strand). Preferred anchor.
  --tss-bed       fallback only when --gtf is not given: BED whose regions are anchored at their 5' end per strand
                  (start for '+', end for '-'); a WARN is logged. Promoter-window BEDs are NOT TSS-centred.
  --threads       accepted for compatibility with the old deepTools step, ignored
  --upstream      bp before the TSS (default 2000)    --downstream  bp after the TSS (default 2000)
  --bin-size      matrix bin, bp (default 50; upstream and downstream must be multiples of it)
  --center-bp     half-width of the "TSS" window, bp (default 250; bins whose centre lies within +-center-bp)
  --flank-bp      width of each outer flank used as background, bp (default 200)
  --strand-col    1-based BED column holding the strand (default: column 4 if it holds +/-, else column 6)
  --no-dedup      testing: keep duplicate TSSs
  --round6        testing: round each bin value to 6 decimals like computeMatrix's matrix file
Same numbers as computeMatrix reference-point --referencePoint TSS --skipZeros --missingDataAsZero (mean per bin),
  followed by the enrichment summary of the mean profile (column means over the kept TSSs):
  tss_signal = maximum of the mean profile inside the centre window; flank_mean = mean of the two outer flanks;
  enrichment = tss_signal / flank_mean; max_over_min = profile maximum / profile minimum (legacy score).
Output: S.tss_enrichment.tsv (header + one row: sample_id n_regions window_up window_down bin_bp tss_signal
  flank_mean enrichment profile_max profile_min max_over_min), S.tss_profile.svg. No matrix file is written.
Exit 2: upstream/downstream not a multiple of the bin size.
Exit 3: missing or empty --bigwig and --gtf/--tss-bed, bigWig cannot be opened, malformed BED/GTF line,
  GTF without usable transcripts, or empty TSS matrix (no TSS covered by signal).
"""
import argparse
import gzip
import os
import re
import sys
import time

import numpy as np

HEADER = ["sample_id", "n_regions", "window_up", "window_down", "bin_bp", "tss_signal", "flank_mean", "enrichment",
          "profile_max", "profile_min", "max_over_min"]
TID_RE = re.compile(r'transcript_id\s+"?([^";]+)"?')
SKIP_PREFIXES = ("#", "track", "browser")


def die(msg, code=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(code)


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_gtf(path):
    """TSS list of (chrom, start0, end, strand) from a GTF: '+' TSS = start0, '-' TSS = end (0-based exclusive end,
    i.e. the 1-based transcript end), as read_bed/build_matrix expect. Transcripts without a 'transcript' feature
    use min exon start / max exon end. Strand from column 7; '.' strands are skipped."""
    tx = {}
    exons = {}
    with open_text(path) as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                die("%s line %d has fewer than 9 tab-separated columns" % (path, ln))
            if f[2] not in ("transcript", "exon"):
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                die("%s line %d has non-integer coordinates" % (path, ln))
            m = TID_RE.search(f[8])
            if not m:
                continue
            key = (f[0], f[6], m.group(1))
            if f[2] == "transcript":
                tx[key] = (start, end)
            else:
                lo, hi = exons.get(key, (start, end))
                exons[key] = (min(lo, start), max(hi, end))
    for key, span in exons.items():
        tx.setdefault(key, span)
    return [(c, s - 1, e, st) for (c, st, _), (s, e) in tx.items() if st in ("+", "-")]


def read_bed(path, strand_col=None):
    """TSS list of (chrom, start, end, strand); strand '+' or '-' from strand_col (default: column 4 if +/-, else
    column 6), otherwise '.'."""
    tss = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(SKIP_PREFIXES):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 3:
                die("%s line %d has fewer than 3 tab-separated columns" % (path, ln))
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                die("%s line %d has non-integer coordinates" % (path, ln))
            col = strand_col or (4 if len(f) >= 4 and f[3] in ("+", "-") else 6)
            strand = f[col - 1] if len(f) >= col else "."
            if strand not in ("+", "-"):
                strand = "."
            tss.append((f[0], start, end, strand))
    return tss


def resolve_chrom(chrom, chroms):
    """Name of chrom in the bigWig chromosome dict, trying the chr-prefix toggle deepTools uses; None if absent."""
    if chrom in chroms:
        return chrom
    if chrom.startswith("chr"):
        cand = chrom[3:]
        if cand == "M":
            cand = "MT"
    else:
        cand = "chr" + chrom
        if cand == "chrMT":
            cand = "chrM"
    return cand if cand in chroms else None


def build_matrix(bw, tss_list, up, down, bin_size, round6=False, dedup=True):
    """Bin matrix (rows = kept TSSs, columns = bins, genome-to-transcript orientation) and the unique TSS count.

    Rows with a zero mean are dropped (skipZeros). TSSs whose chromosome is absent from the bigWig are all-zero rows.
    """
    chroms = bw.chroms()
    nbins = (up + down) // bin_size
    info = []
    seen = set()
    for chrom, start, end, strand in tss_list:
        ref = end if strand == "-" else start
        if strand == "-":
            ws = ref - down
        else:
            ws = ref - up
        res = resolve_chrom(chrom, chroms)
        if dedup:
            key = (res if res is not None else chrom, ref, strand)
            if key in seen:
                continue
            seen.add(key)
        info.append((res, ws, strand))
    n_unique = len(info)
    mat = np.zeros((n_unique, nbins))
    span = nbins * bin_size
    order = sorted(range(n_unique), key=lambda i: (info[i][0] or "", info[i][1]))
    for i in order:
        res, ws, strand = info[i]
        if res is None:
            continue
        a = max(0, ws)
        b = min(int(chroms[res]), ws + span)
        if b <= a:
            continue
        vals = np.asarray(bw.values(res, int(a), int(b)), dtype=float)
        vals = np.where(np.isnan(vals), 0.0, vals)
        window = np.zeros(span)
        window[a - ws:b - ws] = vals
        bins = window.reshape(nbins, bin_size).mean(axis=1)
        if strand == "-":
            bins = bins[::-1]
        if round6:
            bins = np.round(bins, 6)
        mat[i] = bins
    mat = mat[mat.mean(axis=1) != 0]
    return mat, n_unique


def summarise(profile, up, down, bin_size, center, flank):
    """Enrichment summary of the mean profile, as the strings written to the TSV."""
    prof = np.asarray(profile, dtype=float)
    pos = -up + (np.arange(prof.size) + 0.5) * bin_size
    pmax = float(prof.max())
    pmin = float(prof.min())
    centre = np.abs(pos) <= center
    tmax = float(prof[centre].max()) if centre.any() else pmax
    flank_mask = (pos < -up + flank) | (pos > down - flank)
    flankm = float(np.mean(prof[flank_mask])) if flank_mask.any() else 0.0
    enr = "%.4f" % (tmax / flankm) if flankm > 0 else "NA"
    mom = "%.4f" % (pmax / pmin) if pmin > 0 else "NA"
    return {"tss_signal": "%.6g" % tmax, "flank_mean": "%.6g" % flankm, "enrichment": enr,
            "profile_max": "%.6g" % pmax, "profile_min": "%.6g" % pmin, "max_over_min": mom}


def plot_profile(profile, up, bin_size, sample, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["svg.fonttype"] = "none"
    pos = -up + (np.arange(len(profile)) + 0.5) * bin_size
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(pos, profile)
    ax.axvline(0, linestyle="--", color="grey")
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Mean signal")
    ax.set_title("%s: signal around TSS" % sample)
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)


def nonempty(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def main(argv=None):
    t0 = time.time()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample")
    ap.add_argument("--bigwig")
    ap.add_argument("--gtf")
    ap.add_argument("--tss-bed")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--upstream", type=int, default=2000)
    ap.add_argument("--downstream", type=int, default=2000)
    ap.add_argument("--bin-size", type=int, default=50)
    ap.add_argument("--center-bp", type=int, default=250)
    ap.add_argument("--flank-bp", type=int, default=200)
    ap.add_argument("--strand-col", type=int, default=None)
    ap.add_argument("--no-dedup", action="store_true", help="testing: keep duplicate TSSs")
    ap.add_argument("--round6", action="store_true",
                    help="testing: round each bin value to 6 decimals like computeMatrix's matrix file")
    a, extra = ap.parse_known_args(argv)
    if extra:
        sys.stderr.write("WARNING: ignoring unknown arguments: %s\n" % " ".join(extra))

    ref = a.gtf if a.gtf else a.tss_bed
    if not a.sample or not nonempty(a.bigwig) or not nonempty(ref):
        die("missing --sample or empty --bigwig and --gtf/--tss-bed")
    if a.bin_size <= 0 or a.upstream % a.bin_size or a.downstream % a.bin_size or a.upstream + a.downstream <= 0:
        die("--upstream and --downstream must be multiples of --bin-size (%d)" % a.bin_size, code=2)

    import pyBigWig
    bw = pyBigWig.open(a.bigwig)
    if bw is None:
        die("cannot open bigWig %s" % a.bigwig)
    if a.gtf:
        tss = read_gtf(a.gtf)
        if not tss:
            die("no transcripts with strand +/- in %s" % a.gtf)
    else:
        sys.stderr.write("WARN: no --gtf; --tss-bed regions are anchored at their 5' end per strand (start for +, "
                         "end for -), so a promoter-window BED is not centred on the TSS\n")
        tss = read_bed(a.tss_bed, a.strand_col)
    mat, n_unique = build_matrix(bw, tss, a.upstream, a.downstream, a.bin_size, a.round6, not a.no_dedup)
    bw.close()
    if mat.shape[0] == 0:
        die("empty TSS matrix: no TSS overlaps signal")

    profile = mat.mean(axis=0)
    s = summarise(profile, a.upstream, a.downstream, a.bin_size, a.center_bp, a.flank_bp)
    row = [a.sample, str(mat.shape[0]), str(a.upstream), str(a.downstream), str(a.bin_size), s["tss_signal"],
           s["flank_mean"], s["enrichment"], s["profile_max"], s["profile_min"], s["max_over_min"]]
    with open("%s.tss_enrichment.tsv" % a.sample, "w", encoding="utf-8") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.write("\t".join(row) + "\n")
    plot_profile(profile, a.upstream, a.bin_size, a.sample, "%s.tss_profile.svg" % a.sample)
    sys.stderr.write("qc_tss_enrichment: %s: TSSs in %s %d, unique %d, kept %d, %.1f s\n"
                     % (a.sample, "GTF" if a.gtf else "BED", len(tss), n_unique, mat.shape[0], time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
