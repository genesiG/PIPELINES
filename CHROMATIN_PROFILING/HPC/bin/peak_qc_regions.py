#!/usr/bin/env python3
r"""peak_qc_regions.py -- automatic QC regions (loci) for the peak sets of one run, for TRACKS_INI.

Wired into main.nf as process PEAK_QC_REGIONS (modules/local/peak_qc_regions.nf; one call per run), switched by params
run_peak_qc (default true) with the options peak_qc_loci, peak_qc_negative, peak_qc_max_regions, peak_qc_exclude_bed,
peak_qc_shape_marks; its outputs (regions.tsv, regions.detail.tsv, tracks.loci.txt, params.json, log) are published to
results/tracks/peak_qc/ and drawn as track set peak_qc. Implements SPEC_peak_qc_regions.md (sections 2-8). Stdlib only;
pyBigWig is imported lazily and only when bigWigs are needed. Deterministic: no randomness anywhere.

Called by process PEAK_QC_REGIONS (one call per run):
    peak_qc_regions.py --peaksets peaksets.tsv [--classes classes.tsv] [--bigwigs bigwigs.tsv] \
        [--gtf genes.gtf] [--chrom-sizes x.chrom.sizes] [--blacklist x.bed[,y.bed]] [--exclude-bed z.bed[,w.bed]] \
        [--loci "NAME=SPEC[:FLANK],..."] [--negative auto|none|SPEC[:FLANK]] [--config qc_regions.json] \
        [--class-names retained,lost] [--shape-marks K27me3,K27me2] [--nopeak-igg-max-fold F] \
        [--exclude-chrom-regex REGEX] [--spike-prefix PREFIX] [--mito chrM] [--max-regions 15] --outdir DIR

Arguments:
  --peaksets            peaksets.tsv, header id, mark, bed [, role]. role = reference marks the set used for the shape
                        categories of its mark. Default per mark: the base_peakset of its classification, else the first
                        set listed for the mark. Paths are used as given (relative to the working directory).
  --classes             optional classes.tsv, header id, mark, base_peakset, table, contrast [, ref_group]. table = a
                        classification per_region.tsv (plain or .gz) with columns chr, start0, end, name, logFC, logCPM,
                        FDR, class.
  --bigwigs             optional bigwigs.tsv, header sample, mark, bigwig [, condition]. mark IgG = control, every other
                        mark = target. Without it the nopeak_highsignal category and the signal test of the negative are
                        skipped. Rows of the reference condition (ref_group of classes.tsv) are used when both exist.
  --gtf                 GTF for gene coordinates: span of a gene = min start .. max end of all lines whose gene_name is
                        the gene (GTF 1-based start used as-is). Needed for gene loci, the default loci and the negative.
  --chrom-sizes         chrom.sizes. If given it defines the allowed chromosomes and the chromosome ends for shifts.
  --blacklist           comma list of BED files (pipeline key blacklist). A window with > 10% overlap is rejected.
  --exclude-bed         comma list of extra BED masks for nopeak_highsignal (same 10% rule). Default none. Projects
                        that want satellite tiles kept out should pass the general T2T blacklist here.
  --nopeak-igg-max-fold F   optional: reject nopeak tiles whose IgG fold > F. Default off. Needs an IgG row in
                        --bigwigs, else exit 2.
  --exclude-chrom-regex REGEX  contigs matching (re.search) are dropped. Default _|random|Un|alt|EBV|decoy|^chrM$|^MT$
  --spike-prefix PREFIX spike-in contigs starting with PREFIX are dropped. Default none.
  --mito NAME           mitochondrial chromosome, dropped. Default chrM.
  --loci TEXT           user loci, comma list of [NAME=]SPEC[:FLANK]. SPEC = gene symbol, glob (HOXA*; case-insensitive,
                        fnmatch over gene_name), or chr:start-end (used as given, FLANK ignored). Default FLANK 100000.
                        Overrides the config loci. Loci are never filtered and never dropped by --max-regions.
  --negative            auto | none | SPEC[:FLANK] (gene or chr:start-end). Overrides the config negative. Default auto.
  --config FILE         loci and negative in JSON: {"loci": [{"name", "spec", "flank"}], "negative": "none" | "auto" |
                        "SPEC[:FLANK]" | {"name", "spec", "flank"}}. YAML is read only when PyYAML is importable; otherwise
                        exit 2 asking for JSON.
  --class-names A,B     names of the two searched classes, the retained-like one first. Default retained,lost.
  --shape-marks LIST    marks for the shape categories (largest, median_peak, densest). Required when the peaksets have
                        more than 3 marks; at most 3 marks.
  --max-regions N       cap on the number of regions (default 15); the category order decides, user loci always kept.
  --outdir DIR          output directory (created if missing).

Algorithm (SPEC section 4): peak order K = (chrom string, start, end); windows [start, end) with
win_width/centred (500 bp rounding, chromosome-end shift); candidate filters: blacklist <= 10 %, overlap with
accepted windows <= 20 %, at most 200 tries per category (largest: next-best in (-width, K) order; median_peak: next
in n//2, n//2+1, n//2-1, ... order; densest: next-best start); densest = [start_i, start_i + 300 kb) maximising
peak starts (tie: first in K); median = element n//2 of a stable width sort; n < 50 -> no median/densest; class regions
(minority class, rank logCPM desc, FDR asc, width desc, K); cross-mark pairs (50 % overlap of the shorter peak, rank by
min and sum logCPM); nopeak_highsignal and the signal test of the negative only with bigWigs (see class Signal).

Outputs (--outdir):
  regions.tsv           NN_label<TAB>chr:start-end, no header, final order, NN = 01.. over accepted regions.
                        Labels: NN_<locus>, NN_<mark>_largest|median_peak|densest, NN_negative_<GENE|locus>,
                        NN_nopeak_highsignal, NN_<mark>_<class>, NN_<A>_<class>_<B>_<class>.
                        Written start = 1 when the computed start is 0 (make_tracks_ini.py needs start >= 1).
  regions.detail.tsv    columns name (label without NN), category, source, anchor, anchor_name, key_values
                        (k=v;k=v incl. window, width, logCPM, FDR, folds ...), rank_tried, n_rejected, status
                        (ok | skipped:<reason>). Skipped categories are listed with their reason.
  tracks.loci.txt       the accepted windows, comma-joined, one line (the locus_or_gene column of tracks.csv).
  peak_qc_regions.log   warnings and decisions (never aborts for a skipped category).
  params.json           every parameter actually used, the resolved loci and the md5 of every input file.

Exit codes: 0 ok; 2 missing or invalid argument or input file (also: bigWig given but pyBigWig not importable;
--nopeak-igg-max-fold given without an IgG bigWig);
3 no usable peak set (no allowed peak in any set) or a user locus / negative that cannot be resolved in the GTF.
Stdout: none. Stderr: one "ERROR: ..." line on exit 2 / 3.

Deviations from the spec text (decided here, to be reviewed):
  - No YAML in the stdlib: --loci and JSON --config replace the YAML file; YAML only if PyYAML is importable.
  - SPEC may hold '|'-separated alternatives; the default HOXA locus is HOXA*|HOTTIP|HOTAIRM1 (HOTTIP and HOTAIRM1
    included through the alternatives).
  - largest, median_peak and densest retry up to 200 candidates in rank order: largest = (-width, K); median_peak =
    stable width sort, indices n//2, n//2+1, n//2-1, n//2+2, ...; densest = starts ranked by count, skipping anchors
    whose 300 kb window would pass the chromosome end (when --chrom-sizes is given). The first candidate of each is the
    SPEC choice, so the accepted window is unchanged whenever it passes the filters.
  - A user negative (--negative SPEC or config) bypasses the 4.1 filters and, like a user locus, is never dropped by
    --max-regions (prototype behaviour for the filters).
  - Signal folds: target fold of a window = mean over marks of (mean over that mark's bigWigs of f_b); score = mean
    over marks (equal weight per mark) of the mean of log2 f_b; medians of m_b use chromosomes present in every bigWig;
    the allowed chromosomes that some bigWig lacks are listed once in the log.
  - Coverage gates are RELATIVE (changed 2026-10-27 after the full-depth EED run; to be reviewed): rc_b(w) = covered fraction of w in
    bigWig b / median covered fraction of b over ALL allowed tiles (cmed_b). With bamCoverage --skipNonCoveredRegions the
    covered fraction is set by read depth (EED: genome-median tile coverage 0.10-0.46 per target bigWig; no peak-free
    tile reached the former absolute 0.5), so an absolute threshold fails on ordinary and shallow libraries; the
    ratio is depth-free to first order. nopeak tile: min over target bigWigs of rc_b >= NOPEAK_MIN_REL_COVER (every
    target library spreads its reads over at least half its typical breadth: no zero-data tile, no tile whose mean
    comes from a sliver). No-IgG negative: MAX over target bigWigs of rc_b >= NEG_MIN_REL_COVER (some target library
    reads the region at typical breadth; the min would tie the gate to the low signal the negative is chosen for).
    cmed_b <= 0 (a bigWig empty in most tiles) falls back to the mean of the positive tile coverages (WARN). The gate is
    weak for bigWigs that store zero-valued bins: pyBigWig counts a stored zero as covered, so rc_b is ~1 almost
    everywhere for such a file.
  - bigWig rows: IgG rows are always kept (any condition); target rows are kept when their condition is the reference
    group. A mark that loses all its target rows that way is logged as WARN (its signal is then missing).
  - Negative control: the folds and the IgG / coverage gates of the housekeeping and data-driven negatives are computed on
    the gene body +- NEG_FLANK_FREE (the extent of the peak-free test); the whole display window fold is reported as
    win_target_fold. Housekeeping gene accepted when the worst mark fold (max over marks, not the mean) is <=
    NEG_HK_MAX_FOLD; data-driven genes ranked by worst mark fold ascending, then gene name, then start (the recorded
    worst_mark_fold). The negative is a display choice: the gene chosen by these rules is not a representative or
    typical region, and no number derived from it should be reported as typical.
  - Report: the general T2T blacklist is not an input, so only the --blacklist and --exclude-bed fractions are reported.
  - Class shares use the table rows on allowed chromosomes only (spec: contigs dropped before any statistic).
"""
import argparse
import bisect
import fnmatch
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import statistics
import sys
from collections import defaultdict

# ---- constants (SPEC section 8) -------------------------------------------------------------------------------------
ROUND_BP = 500
WIDTH_TIERS = [(10000, 100000), (100000, 200000), (250000, 400000)]
WIDTH_STEP = 100000
DENSEST_WIDTH = 300000
OVERLAP_MAX = 0.20
BLACKLIST_MAX = 0.10
DOMINANT_FRAC = 0.5
CROSS_MIN_OVERLAP = 0.5
CROSS_MIN_PAIRS = 5
NOPEAK_TILE = 250000
NOPEAK_STEP = 50000
NOPEAK_PAD = 10000
NOPEAK_MIN_REL_COVER = 0.5  # min over target bigWigs of covered fraction / that bigWig's genome-median tile coverage
MAX_REGIONS = 15
MAX_TRIES = 200
MIN_PEAKS_SHAPE = 50
MAX_SHAPE_MARKS = 3
NEG_MIN_GENE = 30000
NEG_FLANK_FREE = 20000
NEG_MIN_WIDTH = 200000
NEG_IGG_MIN_FOLD = 0.5
NEG_MIN_REL_COVER = 0.5  # no-IgG negative: max over target bigWigs of the same relative coverage
NEG_HK_MAX_FOLD = 1.5
LOCUS_FLANK = 100000
DEFAULT_EXCLUDE_CHROM = r"_|random|Un|alt|EBV|decoy|^chrM$|^MT$"
DEFAULT_MITO = "chrM"
DEFAULT_CLASS_NAMES = "retained,lost"
DEFAULT_LOCI = [("HOXA", "HOXA*|HOTTIP|HOTAIRM1"), ("HOXD", "HOXD*"), ("CDKN2A", "CDKN2A")]
HOUSEKEEPING = ["GAPDH", "ACTB", "B2M", "PPIA", "RPLP0", "HPRT1", "TBP", "GUSB", "PGK1"]
# Marks of the Polycomb / repressive families (SPEC 4.3 and 4.4a); matched case-insensitively with fnmatch.
POLYCOMB_PATTERNS = ["h3k27me3", "k27me3", "h3k27me2", "k27me2", "ezh2", "suz12", "eed", "cbx*", "ring1b", "h2ak119ub"]
REPRESSIVE_PATTERNS = POLYCOMB_PATTERNS + ["h3k9me2", "h3k9me3", "h3k27me*"]

COORD_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")
COORD_FLANK_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)(?::\d+)?$")
SPEC_FLANK_RE = re.compile(r"^(.+):(\d+)$")
GENE_NAME_RE = re.compile(r'(?:^|[;\s])gene_name "([^"]+)"')
TABLE_COLS = ["chr", "start0", "end", "name", "logFC", "logCPM", "FDR", "class"]


class PeakQCExit(Exception):
    """Raised by die(); main() turns it into the exit code."""

    def __init__(self, code):
        Exception.__init__(self, code)
        self.code = code


def die(msg, code=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    raise PeakQCExit(code)


# ---- window geometry (SPEC 4.0) -------------------------------------------------------------------------------------
def rnd(x, u):
    """Half-up rounding of x to a multiple of u."""
    return int(math.floor(x / u + 0.5)) * u


def win_width(w):
    """Window width for an anchor of width w (SPEC 4.0 tiers)."""
    for mx, ww in WIDTH_TIERS:
        if w <= mx:
            return ww
    return int(math.ceil(2.0 * w / WIDTH_STEP)) * WIDTH_STEP


def shift_into(c, a, b, sizes):
    """Keep the width: shift right if a < 0, shift left if b > chromosome size."""
    w = b - a
    if a < 0:
        a, b = 0, w
    if sizes and c in sizes and b > sizes[c]:
        b = sizes[c]
        a = b - w
    return a, b


def centred(c, s, e, sizes=None, width=None):
    """Window of width win_width(e - s) (or width) centred on the rounded centre of [s, e)."""
    w = width or win_width(e - s)
    ctr = rnd((s + e) / 2.0, ROUND_BP)
    a, b = shift_into(c, ctr - w // 2, ctr - w // 2 + w, sizes)
    return (c, a, b)


def fmt_window(w):
    """chr:start-end for the outputs; start 0 is written as 1 (make_tracks_ini.py needs start >= 1)."""
    c, a, b = w
    return "%s:%d-%d" % (c, a if a > 0 else 1, b)


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return "%.6g" % v
    return str(v)


def sanitize(text):
    """Label-safe name (no whitespace, no colon)."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")


def is_polycomb(mark):
    return any(fnmatch.fnmatchcase(mark.lower(), p) for p in POLYCOMB_PATTERNS)


def is_repressive(mark):
    return any(fnmatch.fnmatchcase(mark.lower(), p) for p in REPRESSIVE_PATTERNS)


# ---- peak selection helpers (SPEC 4.2 - 4.3) ------------------------------------------------------------------------
def largest(peaks):
    """Widest peak; tie -> first in K order."""
    best = None
    for x in peaks:
        if best is None or x[2] - x[1] > best[2] - best[1]:
            best = x
    return best


def median_peak(peaks):
    """Element n // 2 of the stable ascending-by-width sort of the K-sorted peaks (upper median)."""
    s = sorted(peaks, key=lambda x: x[2] - x[1])
    return s[len(s) // 2]


def largest_order(peaks):
    """All peaks by (-width, K): the first is largest(); the rest is the next-best order (SPEC 4.1 rule 2)."""
    return sorted(peaks, key=lambda x: (-(x[2] - x[1]), x[0], x[1], x[2]))


def median_order(peaks):
    """Stable width sort, then n//2, n//2+1, n//2-1, n//2+2, ...; the first is median_peak()."""
    s = sorted(peaks, key=lambda x: x[2] - x[1])
    m = len(s) // 2
    yield s[m]
    for d in range(1, len(s)):
        if m + d < len(s):
            yield s[m + d]
        if m - d >= 0:
            yield s[m - d]


def densest_ranking(peaks, width):
    """All distinct peak starts ranked by (-count, chrom, start); count = peak starts in [start, start + width)."""
    by = defaultdict(list)
    for x in peaks:
        by[x[0]].append(x[1])
    ranking = []
    for c in sorted(by):
        st = by[c]  # ascending, peaks are K-sorted
        prev = None
        for s in st:
            if s == prev:
                continue
            prev = s
            n = bisect.bisect_left(st, s + width) - bisect.bisect_left(st, s)
            ranking.append((n, c, s))
    ranking.sort(key=lambda t: (-t[0], t[1], t[2]))
    return ranking


def frac_overlap(w, chosen):
    """Fraction of window w covered by the union of the chosen windows on the same chromosome."""
    c, a, b = w
    ivs = sorted((x[1], x[2]) for x in chosen if x[0] == c and x[2] > a and x[1] < b)
    cur = a
    ov = 0
    for s, e in ivs:
        s = max(s, cur)
        e = min(e, b)
        if e > s:
            ov += e - s
            cur = e
    return ov / float(b - a)


class Mask:
    """Union of half-open intervals per chromosome for overlap queries."""

    def __init__(self, intervals=()):
        by = defaultdict(list)
        for c, s, e in intervals:
            if e > s:
                by[c].append((s, e))
        self.starts, self.ends = {}, {}
        for c, iv in by.items():
            iv.sort()
            ms, me = [], []
            for s, e in iv:
                if me and s <= me[-1]:
                    me[-1] = max(me[-1], e)
                else:
                    ms.append(s)
                    me.append(e)
            self.starts[c], self.ends[c] = ms, me

    def covered(self, c, a, b):
        ms = self.starts.get(c)
        if not ms or b <= a:
            return 0
        me = self.ends[c]
        i = bisect.bisect_left(ms, b)
        j = bisect.bisect_right(me, a)
        return sum(min(me[k], b) - max(ms[k], a) for k in range(j, i))

    def frac(self, c, a, b):
        return self.covered(c, a, b) / float(b - a)


class Chooser:
    """Accepted windows so far and the candidate filters of SPEC 4.1."""

    def __init__(self, blacklist):
        self.blacklist = blacklist
        self.chosen = []

    def reject_reason(self, w):
        bf = self.blacklist.frac(*w)
        if bf > BLACKLIST_MAX:
            return "blacklist %.3f" % bf
        of = frac_overlap(w, self.chosen)
        if of > OVERLAP_MAX:
            return "overlap %.3f" % of
        return None


def cross_pairs(a_rows, b_rows):
    """Pairs (a in a_rows, b in b_rows) on one chromosome with overlap >= 50 % of the shorter peak, ranked (SPEC 4.7).

    Rows are (chr, start, end, name, logFC, logCPM, FDR, class). Bisect on the sorted starts of b; only the key and the
    two rows are kept per pair."""
    by = defaultdict(list)
    for r in b_rows:
        by[r[0]].append(r)
    index = {}
    for c, lst in by.items():
        lst.sort(key=lambda r: (r[1], r[2]))
        index[c] = ([r[1] for r in lst], lst, max(r[2] - r[1] for r in lst))
    pairs = []
    for a in a_rows:
        if a[0] not in index:
            continue
        starts, lst, maxw = index[a[0]]
        lo = bisect.bisect_right(starts, a[1] - maxw)
        hi = bisect.bisect_left(starts, a[2])
        for k in range(lo, hi):
            b = lst[k]
            ov = min(a[2], b[2]) - max(a[1], b[1])
            if ov <= 0:
                continue
            sh = min(a[2] - a[1], b[2] - b[1])
            if sh <= 0 or ov / float(sh) < CROSS_MIN_OVERLAP:
                continue
            key = (-min(a[5], b[5]), -(a[5] + b[5]), a[0], a[1], a[2], b[1], b[2])
            pairs.append((key, a, b))
    pairs.sort(key=lambda t: t[0])
    return pairs


# ---- readers --------------------------------------------------------------------------------------------------------
def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_tsv(path, required):
    """Tab-separated file with a header line; returns a list of dicts."""
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        missing = [c for c in required if c not in header]
        if missing:
            die("%s: missing column(s) %s" % (path, ", ".join(missing)), 2)
        rows = []
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            rows.append({h: (f[i].strip() if i < len(f) else "") for i, h in enumerate(header)})
    return rows


def read_bed(path):
    """BED3+ as (chrom, start, end, name), sorted by K = (chrom string, start, end)."""
    rows = []
    with open(path, "r") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\r\n").split("\t")
            try:
                rows.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else ""))
            except (ValueError, IndexError):
                die("%s line %d: BED start/end not integers" % (path, ln), 2)
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    return rows


def read_table(path):
    """classification per_region table as (chr, start0, end, name, logFC, logCPM, FDR, class), K-sorted."""
    rows = []
    with open_text(path) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        missing = [c for c in TABLE_COLS if c not in header]
        if missing:
            die("%s: missing column(s) %s" % (path, ", ".join(missing)), 2)
        ix = {c: header.index(c) for c in TABLE_COLS}
        for ln, line in enumerate(fh, 2):
            if not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            try:
                rows.append((f[ix["chr"]], int(f[ix["start0"]]), int(f[ix["end"]]), f[ix["name"]],
                             float(f[ix["logFC"]]), float(f[ix["logCPM"]]), float(f[ix["FDR"]]), f[ix["class"]]))
            except (ValueError, IndexError):
                die("%s line %d: bad value in %s" % (path, ln, ", ".join(TABLE_COLS)), 2)
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    return rows


def read_sizes(path):
    sizes = {}
    with open(path, "r") as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                sizes[f[0]] = int(f[1])
    return sizes


def read_gene_spans(path):
    """gene_name -> [set of chromosomes, min start, max end] over all GTF lines carrying gene_name."""
    spans = {}
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("#") or 'gene_name "' not in line:
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                continue
            m = GENE_NAME_RE.search(f[8])
            if not m:
                continue
            try:
                s, e = int(f[3]), int(f[4])
            except ValueError:
                die("%s: bad start/end in GTF line" % path, 2)
            v = spans.get(m.group(1))
            if v is None:
                spans[m.group(1)] = [set([f[0]]), s, e]
            else:
                v[0].add(f[0])
                v[1] = min(v[1], s)
                v[2] = max(v[2], e)
    return spans


def lookup_gene(spans, pattern):
    """(chrom, start, end) of the span over all gene_name values matching the symbol or glob (case-insensitive)."""
    pats = [p.upper() for p in pattern.split("|") if p]  # "A*|B" = alternatives
    hits = [v for name, v in spans.items() if any(fnmatch.fnmatchcase(name.upper(), p) for p in pats)]
    if not hits:
        return None
    chroms = set()
    for h in hits:
        chroms |= h[0]
    if len(chroms) != 1:
        raise LocusError("%s spans several chromosomes (%s)" % (pattern, ",".join(sorted(chroms))))
    return (chroms.pop(), min(h[1] for h in hits), max(h[2] for h in hits))


class LocusError(Exception):
    pass


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---- loci and negative specs (SPEC 2, 3.4) --------------------------------------------------------------------------
def parse_spec(text, default_flank):
    """'chr:s-e[:x]' -> coord; 'SPEC:FLANK' -> spec + flank; 'SPEC' -> spec + default flank."""
    text = text.strip()
    if not text:
        die("empty locus specification", 2)
    m = COORD_FLANK_RE.match(text)
    if m:
        c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        if e <= s:
            die("locus %s: end must be greater than start" % text, 2)
        return {"spec": "%s:%d-%d" % (c, s, e), "coord": (c, s, e), "flank": None}
    m = SPEC_FLANK_RE.match(text)
    if m:
        return {"spec": m.group(1), "coord": None, "flank": int(m.group(2))}
    return {"spec": text, "coord": None, "flank": default_flank}


def make_locus(item, default_flank=LOCUS_FLANK, name=None):
    """One locus from 'NAME=SPEC[:FLANK]' (name optional)."""
    item = item.strip()
    if "=" in item and name is None:
        name, item = item.split("=", 1)
        name = name.strip()
    p = parse_spec(item, default_flank)
    if not name:
        name = sanitize(p["spec"])
    return {"name": sanitize(name), "spec": p["spec"], "coord": p["coord"], "flank": p["flank"]}


def parse_loci_text(text):
    loci = []
    for item in text.split(","):
        if item.strip():
            loci.append(make_locus(item))
    return loci


def load_config(path):
    with open(path, "r") as fh:
        text = fh.read()
    try:
        data = json.loads(text)
    except ValueError:
        try:
            import yaml
        except ImportError:
            die("config %s is not JSON and PyYAML is not importable; write the config as JSON" % path, 2)
        try:
            data = yaml.safe_load(text)
        except Exception as ex:  # yaml errors are not a common base class
            die("config %s cannot be parsed: %s" % (path, ex), 2)
    if data is None:
        data = {}
    if not isinstance(data, dict):
        die("config %s must be a mapping with keys loci and negative" % path, 2)
    return data


def loci_from_config(entries):
    loci = []
    for e in entries or []:
        if isinstance(e, str):
            loci.append(make_locus(e))
        elif isinstance(e, dict) and e.get("spec"):
            flank = e.get("flank")
            loci.append(make_locus(str(e["spec"]) if e.get("name") is None else
                                   "%s=%s" % (e["name"], e["spec"]),
                                   LOCUS_FLANK if flank is None else int(flank)))
        else:
            die("config loci entry needs a spec: %r" % (e,), 2)
    return loci


def negative_from_text(text):
    v = text.strip()
    if v.lower() in ("", "auto"):
        return "auto"
    if v.lower() == "none":
        return "none"
    return make_locus(v)


# ---- signal (SPEC 3.3, 4.5) -----------------------------------------------------------------------------------------
def window_stats(bw, chrom, a, b):
    """(full-width mean, covered fraction) of [a, b) on a pyBigWig handle; bases without data count as 0.

    pyBigWig 'mean' with exact=True averages over the covered bases only, so the full-width mean is
    mean * coverage (coverage = covered fraction of the requested range); the clipped part of the window to the
    chromosome counts as uncovered."""
    size = bw.chroms(chrom)
    if not size or b <= a:
        return 0.0, 0.0
    a2, b2 = max(a, 0), min(b, size)
    if b2 <= a2:
        return 0.0, 0.0
    mean = bw.stats(chrom, a2, b2, type="mean", nBins=1, exact=True)[0]
    cov = bw.stats(chrom, a2, b2, type="coverage", nBins=1, exact=True)[0]
    if mean is None or cov is None or mean != mean:
        return 0.0, 0.0
    frac_in = (b2 - a2) / float(b - a)
    return mean * cov * frac_in, cov * frac_in


def window_mean(bw, chrom, a, b):
    """Full-width mean of [a, b) (missing bins = 0), SPEC 4.5.3."""
    return window_stats(bw, chrom, a, b)[0]


class Signal:
    """Target and IgG bigWigs: per-window folds f_b and the tile medians med_b (SPEC 4.5)."""

    def __init__(self, rows, sizes, allowed, log):
        self.rows = rows
        self.sizes = sizes
        self.allowed = allowed
        self.log = log
        self.bws = []
        self.tiles = None
        self.vals = None
        self.med = None
        self.cmed = None
        self.target = [i for i, r in enumerate(rows) if r["mark"].lower() != "igg"]
        self.igg = [i for i, r in enumerate(rows) if r["mark"].lower() == "igg"]

    def open(self):
        if self.bws:
            return
        try:
            import pyBigWig
        except ImportError:
            die("--bigwigs given but pyBigWig is not importable (use the cutandrun-py container)", 2)
        for r in self.rows:
            try:
                bw = pyBigWig.open(r["bigwig"])
                ok = bw is not None and bw.isBigWig()
            except RuntimeError:
                ok = False
            if not ok:
                die("not a bigWig: %s" % r["bigwig"], 2)
            self.bws.append(bw)

    def window(self, c, a, b):
        self.open()
        return [window_stats(bw, c, a, b) for bw in self.bws]

    def tile_table(self):
        """All allowed tiles (250 kb, step 50 kb, inside the chromosome) with their (mean, coverage) per bigWig, the
        median of the means (med) and the median of the covered fractions (cmed) per bigWig over ALL tiles."""
        if self.tiles is None:
            self.open()
            names = self.sizes.keys() if self.sizes is not None else self.bws[0].chroms().keys()
            chroms = sorted(c for c in names if self.allowed(c) and all(bw.chroms(c) for bw in self.bws))
            dropped = sorted(c for c in names if self.allowed(c) and not all(bw.chroms(c) for bw in self.bws))
            if dropped:
                self.log("signal: allowed chromosomes dropped because a bigWig lacks them: %s" % ",".join(dropped))
            size_of = (lambda c: self.sizes[c]) if self.sizes is not None else (lambda c: self.bws[0].chroms(c))
            tiles, vals = [], []
            for c in chroms:
                for a in range(0, size_of(c) - NOPEAK_TILE + 1, NOPEAK_STEP):
                    tiles.append((c, a, a + NOPEAK_TILE))
                    vals.append(self.window(c, a, a + NOPEAK_TILE))
            self.tiles, self.vals = tiles, vals
            self.med = []
            for i, r in enumerate(self.rows):
                xs = [v[i][0] for v in vals]
                med = statistics.median(xs) if xs else 0.0
                if med <= 0:
                    pos = [x for x in xs if x > 0]
                    med = statistics.mean(pos) if pos else 1.0
                    self.log("WARN median of bigWig %s is 0; fold uses %.6g" % (r["bigwig"], med))
                self.med.append(med)
            self.cmed = []
            for i, r in enumerate(self.rows):
                cs = [v[i][1] for v in vals]
                cmed = statistics.median(cs) if cs else 0.0
                if cmed <= 0:
                    pos = [x for x in cs if x > 0]
                    cmed = statistics.mean(pos) if pos else 1.0
                    self.log("WARN median tile coverage of bigWig %s is 0; relative coverage uses %.6g"
                             % (r["bigwig"], cmed))
                self.cmed.append(cmed)
            self.log("signal: %d tiles on %d chromosomes" % (len(tiles), len(chroms)))
            self.log("signal: median tile mean / covered fraction per bigWig: %s" % ", ".join(
                "%s %.4g/%.4g" % (r["sample"], self.med[i], self.cmed[i]) for i, r in enumerate(self.rows)))
        return self.tiles, self.vals

    def summary(self, vals):
        """Folds and relative coverages of one window from its (mean, coverage) list (SPEC 4.5.2 - 4.5.4).
        cov = min absolute covered fraction over the target bigWigs (reported only); rel_cov_min / rel_cov_max = min /
        max over the target bigWigs of covered fraction / cmed_b (the gates)."""
        self.tile_table()
        folds = [(vals[i][0] + 0.05 * self.med[i]) / (1.05 * self.med[i]) for i in range(len(vals))]
        rel = [vals[i][1] / self.cmed[i] for i in self.target]
        by_mark = defaultdict(list)
        for i in self.target:
            by_mark[self.rows[i]["mark"]].append(i)
        mark_fold = {}
        mark_score = []
        for m, ids in by_mark.items():
            mark_fold[m] = sum(folds[i] for i in ids) / len(ids)
            mark_score.append(sum(math.log(folds[i], 2) for i in ids) / len(ids))
        igg = [folds[i] for i in self.igg]
        return {
            "mark_fold": mark_fold,
            "target_fold": (sum(mark_fold.values()) / len(mark_fold)) if mark_fold else None,
            "igg_fold": (sum(igg) / len(igg)) if igg else None,
            "cov": min(vals[i][1] for i in self.target) if self.target else None,
            "rel_cov_min": min(rel) if rel else None,
            "rel_cov_max": max(rel) if rel else None,
            "score": (sum(mark_score) / len(mark_score)) if mark_score else None,
        }

    def kv(self, sm):
        kv = {}
        for m, f in sm["mark_fold"].items():
            kv["fold_" + sanitize(m)] = f
        kv["target_fold"] = sm["target_fold"]
        kv["igg_fold"] = sm["igg_fold"]
        kv["cov_min"] = sm["cov"]
        kv["rel_cov_min"] = sm["rel_cov_min"]
        kv["rel_cov_max"] = sm["rel_cov_max"]
        kv["score"] = sm["score"]
        return kv


# ---- category builders ----------------------------------------------------------------------------------------------
def skipped(name, category, source, reason, kv=None):
    return {"name": name, "category": category, "source": source, "anchor": "-", "anchor_name": "-",
            "kv": dict(kv or {}), "rank": 0, "rej": 0, "status": "skipped:" + reason, "window": None}


def forced(name, category, source, w, info, chooser, user=False):
    """User locus or user negative: accepted without the 4.1 filters, before ranked categories. user = never dropped
    by --max-regions (user loci only)."""
    chooser.chosen.append(w)
    kv = dict(info.get("kv", {}))
    kv["window"] = fmt_window(w)
    kv["window_width"] = w[2] - w[1]
    return {"name": info.get("name", name), "category": category, "source": source,
            "anchor": info.get("anchor", "-"), "anchor_name": info.get("anchor_name", "-"), "kv": kv,
            "rank": 1, "rej": 0, "status": "ok", "window": w, "user": user}


def place(name, category, source, cands, chooser, kv0=None):
    """First candidate (window, info) in rank order that passes the 4.1 filters (at most MAX_TRIES)."""
    n_rej = 0
    last = None
    tried = 0
    for k, (w, info) in enumerate(itertools.islice(cands, MAX_TRIES), 1):
        tried = k
        why = chooser.reject_reason(w)
        if why is None:
            chooser.chosen.append(w)
            kv = dict(kv0 or {})
            kv.update(info.get("kv", {}))
            kv["window"] = fmt_window(w)
            kv["window_width"] = w[2] - w[1]
            return {"name": info.get("name", name), "category": category, "source": source,
                    "anchor": info.get("anchor", "-"), "anchor_name": info.get("anchor_name", "-"), "kv": kv,
                    "rank": k, "rej": n_rej, "status": "ok", "window": w}
        n_rej += 1
        last = why
    reason = "no candidate" if tried == 0 else "no candidate passed filters (%d tried; last %s)" % (tried, last)
    return {"name": name, "category": category, "source": source, "anchor": "-", "anchor_name": "-",
            "kv": dict(kv0 or {}), "rank": tried, "rej": n_rej, "status": "skipped:" + reason, "window": None}


def peak_info(x, extra=None):
    kv = {"width": x[2] - x[1]}
    kv.update(extra or {})
    return {"anchor": "%s:%d-%d" % (x[0], x[1], x[2]), "anchor_name": x[3], "kv": kv}


def class_info(r):
    return {"anchor": "%s:%d-%d" % (r[0], r[1], r[2]), "anchor_name": r[3],
            "kv": {"logFC": r[4], "logCPM": r[5], "FDR": r[6], "width": r[2] - r[1]}}


def nopeak_category(sig, pad_mask, bl, exclude, igg_max, chooser, log):
    """nopeak_highsignal: highest-scoring eligible 250 kb tile (SPEC 4.5)."""
    name, category = "nopeak_highsignal", "nopeak"
    source = "bigwigs:" + ",".join(r["sample"] for r in sig.rows)
    if not sig.target:
        return skipped(name, category, source, "no target bigWig")
    tiles, vals = sig.tile_table()
    if not tiles:
        return skipped(name, category, source, "no tile on the allowed chromosomes")
    cands = []
    n_free = n_mask = n_cov = 0
    for (c, a, b), v in zip(tiles, vals):
        if pad_mask.covered(c, a, b) > 0:
            continue
        n_free += 1
        if bl.frac(c, a, b) > BLACKLIST_MAX:
            continue
        if exclude is not None and exclude.frac(c, a, b) > BLACKLIST_MAX:
            continue
        n_mask += 1
        sm = sig.summary(v)
        if sm["rel_cov_min"] < NOPEAK_MIN_REL_COVER:
            continue
        n_cov += 1
        if igg_max is not None and sm["igg_fold"] > igg_max:
            continue  # igg_max is only set when an IgG bigWig exists (checked once in run)
        cands.append((-sm["score"], c, a, b, sm))
    counts = "%d tiles; %d peak-free; %d after blacklist/exclude; %d with relative coverage >= %g; %d after IgG cap" % (
        len(tiles), n_free, n_mask, n_cov, NOPEAK_MIN_REL_COVER, len(cands))
    log("nopeak_highsignal: " + counts)
    if not cands:
        return skipped(name, category, source, "no eligible tile (%s)" % counts)
    cands.sort(key=lambda t: (t[0], t[1], t[2], t[3]))

    def gen():
        for negscore, c, a, b, sm in cands:
            kv = sig.kv(sm)
            kv["n_eligible"] = len(cands)
            kv["bl_frac"] = bl.frac(c, a, b)
            if exclude is not None:
                kv["exclude_frac"] = exclude.frac(c, a, b)
            yield (c, a, b), {"anchor": "%s:%d-%d" % (c, a, b), "anchor_name": "-", "kv": kv}

    return place(name, category, source, gen(), chooser)


def shape_category(mark, src_id, P, n_all, sizes, chooser):
    recs = []
    if not P:
        for suffix in ("largest", "median_peak", "densest"):
            recs.append(skipped("%s_%s" % (mark, suffix), "shape", src_id, "no allowed peak"))
        return recs
    def peak_cands(order):
        for x in order:
            yield centred(x[0], x[1], x[2], sizes), peak_info(x, {"n": n_all})

    recs.append(place("%s_largest" % mark, "shape", src_id, peak_cands(largest_order(P)), chooser))
    if n_all < MIN_PEAKS_SHAPE:
        for suffix in ("median_peak", "densest"):
            recs.append(skipped("%s_%s" % (mark, suffix), "shape", src_id, "n=%d<%d" % (n_all, MIN_PEAKS_SHAPE)))
        return recs
    recs.append(place("%s_median_peak" % mark, "shape", src_id, peak_cands(median_order(P)), chooser))
    ranking = densest_ranking(P, DENSEST_WIDTH)

    def gen():
        for n_starts, c, s in ranking:
            if sizes and c in sizes and s + DENSEST_WIDTH > sizes[c]:
                continue  # the 300 kb window would pass the chromosome end
            yield ((c, s, s + DENSEST_WIDTH),
                   {"anchor": "%s:%d-%d" % (c, s, s + DENSEST_WIDTH), "anchor_name": "-",
                    "kv": {"n_starts": n_starts, "n": n_all}})

    recs.append(place("%s_densest" % mark, "shape", src_id, gen(), chooser))
    return recs


def class_category(cl, ret, lost, sizes, chooser):
    recs = []
    rows = cl["rows"]
    tot = len(rows)
    for cname in (ret, lost):
        name = "%s_%s" % (cl["mark"], cname)
        members = [r for r in rows if r[7] == cname]
        n_c = len(members)
        share = n_c / float(tot) if tot else 0.0
        kv0 = {"n_class": n_c, "n_total": tot, "share": share}
        if n_c == 0:
            recs.append(skipped(name, "class", cl["id"], "empty class (n=0)", kv0))
            continue
        if share > DOMINANT_FRAC:
            recs.append(skipped(name, "class", cl["id"], "dominant class share=%.3f" % share, kv0))
            continue
        members.sort(key=lambda r: (-r[5], r[6], -(r[2] - r[1]), r[0], r[1], r[2]))

        def gen(members=members):
            for r in members:
                yield centred(r[0], r[1], r[2], sizes), class_info(r)

        recs.append(place(name, "class", cl["id"], gen(), chooser, kv0))
    return recs


def cross_category(A, B, a_cls, b_cls, sizes, chooser):
    """One orientation: A rows of class a_cls against B rows of class b_cls (SPEC 4.7)."""
    recs = []
    name = "%s_%s_%s_%s" % (A["mark"], a_cls, B["mark"], b_cls)
    a_rows = [r for r in A["rows"] if r[7] == a_cls]
    b_rows = [r for r in B["rows"] if r[7] == b_cls]
    pairs = cross_pairs(a_rows, b_rows)
    src = "%s/%s" % (A["id"], B["id"])
    kv0 = {"n_pairs": len(pairs)}
    if len(pairs) < CROSS_MIN_PAIRS:
        return [skipped(name, "cross", src, "pairs=%d<%d" % (len(pairs), CROSS_MIN_PAIRS), kv0)]

    def gen():
        for _, a, b in pairs:
            c = a[0]
            s, e = min(a[1], b[1]), max(a[2], b[2])
            yield centred(c, s, e, sizes), {
                "anchor": "%s:%d-%d" % (c, s, e), "anchor_name": "%s|%s" % (a[3], b[3]),
                "kv": {"cpm_%s" % A["mark"]: a[5], "cpm_%s" % B["mark"]: b[5], "fdr_%s" % A["mark"]: a[6],
                       "fdr_%s" % B["mark"]: b[6], "n_pairs": len(pairs)}}

    recs.append(place(name, "cross", src, gen(), chooser, kv0))
    return recs


# ---- run ------------------------------------------------------------------------------------------------------------
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Automatic peak-QC regions (see the module docstring).")
    ap.add_argument("--peaksets", required=True)
    ap.add_argument("--classes", default=None)
    ap.add_argument("--bigwigs", default=None)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--chrom-sizes", default=None)
    ap.add_argument("--blacklist", default="")
    ap.add_argument("--exclude-bed", default="")
    ap.add_argument("--nopeak-igg-max-fold", type=float, default=None)
    ap.add_argument("--exclude-chrom-regex", default=DEFAULT_EXCLUDE_CHROM)
    ap.add_argument("--spike-prefix", default="")
    ap.add_argument("--mito", default=DEFAULT_MITO)
    ap.add_argument("--loci", default=None)
    ap.add_argument("--negative", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--class-names", default=DEFAULT_CLASS_NAMES)
    ap.add_argument("--shape-marks", default=None)
    ap.add_argument("--max-regions", type=int, default=MAX_REGIONS)
    ap.add_argument("--outdir", required=True)
    return ap.parse_args(argv)


def split_list(text):
    return [t.strip() for t in (text or "").split(",") if t.strip()]


def run(args, log, inputs):
    """Fill inputs (list of input paths for the md5s) and return the category records and metadata."""
    for label, path in [("--peaksets", args.peaksets), ("--classes", args.classes), ("--bigwigs", args.bigwigs),
                        ("--gtf", args.gtf), ("--chrom-sizes", args.chrom_sizes), ("--config", args.config)]:
        if path and not os.path.isfile(path):
            die("%s file not found: %s" % (label, path), 2)
    blacklist_files = split_list(args.blacklist)
    exclude_files = split_list(args.exclude_bed)
    inputs.extend(blacklist_files + exclude_files)
    for path in blacklist_files + exclude_files:
        if not os.path.isfile(path):
            die("BED file not found: %s" % path, 2)
    if args.max_regions < 1:
        die("--max-regions must be >= 1", 2)
    class_names = split_list(args.class_names)
    if len(class_names) != 2:
        die("--class-names needs two names (retained-like,lost-like)", 2)
    RET, LOST = class_names
    try:
        re.compile(args.exclude_chrom_regex)
    except re.error as ex:
        die("--exclude-chrom-regex invalid: %s" % ex, 2)
    cfg = load_config(args.config) if args.config else {}

    sizes = read_sizes(args.chrom_sizes) if args.chrom_sizes else None
    exclude_rx = re.compile(args.exclude_chrom_regex)

    def allowed(c):
        if sizes is not None and c not in sizes:
            return False
        if args.mito and c == args.mito:
            return False
        if args.spike_prefix and c.startswith(args.spike_prefix):
            return False
        return exclude_rx.search(c) is None

    spans_box = {}

    def spans():
        if "s" not in spans_box:
            if not args.gtf:
                die("gene locus, default loci or negative need --gtf", 2)
            spans_box["s"] = read_gene_spans(args.gtf)
        return spans_box["s"]

    # peak sets
    ps_rows = read_tsv(args.peaksets, ["id", "mark", "bed"])
    peaksets = {}
    mark_order = []
    for r in ps_rows:
        if not os.path.isfile(r["bed"]):
            die("BED file not found: %s" % r["bed"], 2)
        if r["mark"] not in mark_order:
            mark_order.append(r["mark"])
        inputs.append(r["bed"])
        peaks = [x for x in read_bed(r["bed"]) if allowed(x[0])]
        peaksets[r["id"]] = {"id": r["id"], "mark": r["mark"], "role": r.get("role", ""), "peaks": peaks,
                             "bed": r["bed"], "n_raw": None}
        log("peak set %s (%s): %d allowed peaks" % (r["id"], r["mark"], len(peaks)))
    if not any(p["peaks"] for p in peaksets.values()):
        die("no usable peak set (no allowed peak in %s)" % args.peaksets, 3)

    # classifications
    classes = []
    if args.classes:
        for r in read_tsv(args.classes, ["id", "mark", "base_peakset", "table", "contrast"]):
            if r["base_peakset"] not in peaksets:
                die("classes.tsv %s: base_peakset %s is not in peaksets.tsv" % (r["id"], r["base_peakset"]), 2)
            if not os.path.isfile(r["table"]):
                die("classification table not found: %s" % r["table"], 2)
            inputs.append(r["table"])
            rows = [x for x in read_table(r["table"]) if allowed(x[0])]
            classes.append({"id": r["id"], "mark": r["mark"], "base": r["base_peakset"], "table": r["table"],
                            "contrast": r["contrast"], "ref_group": r.get("ref_group", ""), "rows": rows})
            log("classification %s (%s, contrast %s): %d allowed rows" % (r["id"], r["mark"], r["contrast"], len(rows)))

    # reference peak set per mark
    ref = {}
    for m in mark_order:
        in_mark = [p for p in peaksets.values() if p["mark"] == m]
        chosen = [p for p in in_mark if p["role"] == "reference"]
        if not chosen:
            chosen = [p for p in peaksets.values() if p["id"] in [c["base"] for c in classes if c["mark"] == m]]
        if not chosen:
            chosen = in_mark
        if len(chosen) > 1:
            log("WARN mark %s has %d reference sets; the first is used" % (m, len(chosen)))
        ref[m] = chosen[0]

    # shape marks
    if args.shape_marks:
        shape = split_list(args.shape_marks)
        for m in shape:
            if m not in mark_order:
                die("--shape-marks %s is not a mark of peaksets.tsv" % m, 2)
    else:
        shape = mark_order
        if len(shape) > MAX_SHAPE_MARKS:
            die("%d marks: give --shape-marks (at most %d)" % (len(shape), MAX_SHAPE_MARKS), 2)
    if len(shape) > MAX_SHAPE_MARKS:
        die("at most %d shape marks" % MAX_SHAPE_MARKS, 2)

    # masks
    all_peaks = [x for p in peaksets.values() for x in p["peaks"]]
    peak_mask = Mask((x[0], x[1], x[2]) for x in all_peaks)
    pad_mask = Mask((x[0], max(0, x[1] - NOPEAK_PAD), x[2] + NOPEAK_PAD) for x in all_peaks)
    blacklist = Mask(iv for f in blacklist_files for iv in ((x[0], x[1], x[2]) for x in read_bed(f)))
    exclude = (Mask(iv for f in exclude_files for iv in ((x[0], x[1], x[2]) for x in read_bed(f)))
               if exclude_files else None)
    chooser = Chooser(blacklist)
    recs = []

    # bigWigs
    bw_rows = []
    if args.bigwigs:
        bw_rows = read_tsv(args.bigwigs, ["sample", "mark", "bigwig"])
        for r in bw_rows:
            if not os.path.isfile(r["bigwig"]):
                die("bigWig not found: %s" % r["bigwig"], 2)
            inputs.append(r["bigwig"])
        ref_groups = set(c["ref_group"] for c in classes if c["ref_group"])
        if ref_groups:
            # IgG rows are kept whatever their condition: one control serves all conditions
            kept = [r for r in bw_rows if r["mark"].lower() == "igg" or r.get("condition", "") in ref_groups]
            if any(r["mark"].lower() != "igg" for r in kept):
                had = set(r["mark"] for r in bw_rows if r["mark"].lower() != "igg")
                now = set(r["mark"] for r in kept if r["mark"].lower() != "igg")
                for m in sorted(had - now):
                    if m in mark_order:
                        log("WARN mark %s lost all its target bigWigs (no row of reference condition %s); "
                            "it has no signal" % (m, ",".join(sorted(ref_groups))))
                bw_rows = kept
            else:
                log("WARN no bigWig of the reference condition (%s); all bigWigs used" % ",".join(sorted(ref_groups)))
    if args.nopeak_igg_max_fold is not None and not any(r["mark"].lower() == "igg" for r in bw_rows):
        die("--nopeak-igg-max-fold needs an IgG bigWig (mark IgG) in --bigwigs", 2)
    sig_box = {}

    def signal():
        if "sig" not in sig_box:
            sig_box["sig"] = Signal(bw_rows, sizes, allowed, log)
        return sig_box["sig"]

    # 1. user loci or default marker loci
    loci_source = "none"
    loci = []
    if args.loci is not None:
        loci = parse_loci_text(args.loci)
        loci_source = "cli"
    elif cfg.get("loci"):
        loci = loci_from_config(cfg.get("loci"))
        loci_source = "config"
    if loci:
        for loc in loci:
            try:
                if loc["coord"]:
                    w = loc["coord"][0], loc["coord"][1], loc["coord"][2]
                else:
                    g = lookup_gene(spans(), loc["spec"])
                    if g is None:
                        die("locus %s (%s): gene not found in the GTF" % (loc["name"], loc["spec"]), 3)
                    w = (g[0], max(1, g[1] - loc["flank"]), g[2] + loc["flank"])
            except LocusError as ex:
                die("locus %s: %s" % (loc["name"], ex), 3)
            recs.append(forced(loc["name"], "user_locus", "loci", w,
                               {"anchor": "%s:%d-%d" % w, "anchor_name": loc["spec"],
                                "kv": {"spec": loc["spec"], "flank": loc["flank"] if loc["flank"] is not None else "-"}},
                               chooser, user=True))
    elif args.gtf and any(is_polycomb(m) for m in mark_order):
        loci_source = "default"
        kept = 0
        for name, pattern in DEFAULT_LOCI:
            g = None
            try:
                g = lookup_gene(spans(), pattern)
            except LocusError as ex:
                log("WARN default locus %s: %s" % (name, ex))
            if g is None:
                log("default locus %s: not found in the GTF; not used" % name)
                recs.append(skipped(name, "user_locus", "default", "gene not found in GTF"))
                continue
            w = (g[0], max(1, g[1] - LOCUS_FLANK), g[2] + LOCUS_FLANK)
            poly_hit = peak_mask_poly_overlap(w, ref, mark_order)
            if not poly_hit or kept >= MAX_SHAPE_MARKS:
                why = "not a target in this data" if not poly_hit else "max 3 default loci"
                log("default locus %s %s: %s" % (name, fmt_window(w), why))
                recs.append(skipped(name, "user_locus", "default", why))
                continue
            kept += 1
            recs.append(forced(name, "user_locus", "default", w,
                               {"anchor": "%s:%d-%d" % w, "anchor_name": pattern,
                                "kv": {"spec": pattern, "flank": LOCUS_FLANK}}, chooser, user=True))
    else:
        log("no user loci and no default loci (no --gtf or no Polycomb-type mark)")

    # 2. shape categories per mark
    for m in shape:
        p = ref[m]
        recs.extend(shape_category(m, p["id"], p["peaks"], len(p["peaks"]), sizes, chooser))

    # 3. negative control
    neg_text = args.negative
    neg_source = "cli"
    if neg_text is None:
        neg_text = cfg.get("negative")
        neg_source = "config"
    if neg_text is None:
        neg_text, neg_source = "auto", "default"
    if isinstance(neg_text, dict):
        neg = make_locus(str(neg_text.get("spec")), LOCUS_FLANK if neg_text.get("flank") is None
                         else int(neg_text["flank"]), neg_text.get("name"))
    else:
        neg = negative_from_text(str(neg_text))
    recs.extend(negative_category(neg, neg_source, spans, sizes, chooser, peak_mask, blacklist, mark_order,
                                  bw_rows, signal, allowed, log, bool(args.gtf)))

    # 4. nopeak_highsignal
    if bw_rows:
        recs.append(nopeak_category(signal(), pad_mask, blacklist, exclude, args.nopeak_igg_max_fold, chooser, log))
    else:
        recs.append(skipped("nopeak_highsignal", "nopeak", "bigwigs", "no --bigwigs"))

    # 5. classifications: minority classes
    for cl in classes:
        recs.extend(class_category(cl, RET, LOST, sizes, chooser))

    # 6. cross-mark discordant regions
    for i, A in enumerate(classes):
        for B in classes[i + 1:]:
            if A["contrast"] != B["contrast"] or A["mark"] == B["mark"]:
                continue
            if mark_order.index(A["mark"]) > mark_order.index(B["mark"]):
                A, B = B, A
            for a_cls, b_cls in ((LOST, RET), (RET, LOST)):
                recs.extend(cross_category(A, B, a_cls, b_cls, sizes, chooser))
    return recs, loci_source, loci, neg, neg_source, ref, shape


def peak_mask_poly_overlap(w, ref, mark_order):
    """True when a peak of the reference set of a Polycomb-type mark overlaps w (SPEC 4.3)."""
    for m in mark_order:
        if is_polycomb(m):
            for x in ref[m]["peaks"]:
                if x[0] == w[0] and x[2] > w[1] and x[1] < w[2]:
                    return True
    return False


def negative_category(neg, neg_source, spans, sizes, chooser, peak_mask, blacklist, mark_order, bw_rows, signal,
                      allowed, log, has_gtf):
    """negative control (SPEC 4.4). neg: 'none' | 'auto' | locus dict (user negative, bypasses the 4.1 filters).
    The auto negative is a display choice (a peak-free gene body, chosen by the rules below), not a sample."""
    name_default = "negative"
    if neg == "none":
        return [skipped(name_default, "negative", neg_source, "negative none")]
    if isinstance(neg, dict):
        try:
            g = neg["coord"] if neg["coord"] else lookup_gene(spans(), neg["spec"])
        except LocusError as ex:
            die("negative %s: %s" % (neg["spec"], ex), 3)
        if g is None:
            die("negative %s: gene not found in the GTF" % neg["spec"], 3)
        c, s, e = g
        flank = neg["flank"] or 0
        a = max(1, s - flank) if not neg["coord"] else s
        b = e + flank if not neg["coord"] else e
        if b - a < NEG_MIN_WIDTH:
            w = centred(c, a, b, sizes, width=NEG_MIN_WIDTH)
        else:
            w = (c, a, b)
        kv = {"spec": neg["spec"], "flank": neg["flank"], "gene_span": "%s:%d-%d" % (c, s, e)}
        if bw_rows:
            sig = signal()
            kv.update(sig.kv(sig.summary(sig.window(*w))))
        return [forced("negative_" + neg["name"], "negative", neg_source, w,
                       {"name": "negative_" + neg["name"], "anchor": "%s:%d-%d" % w,
                        "anchor_name": neg["spec"], "kv": kv}, chooser, user=True)]

    # auto
    if not has_gtf:
        return [skipped(name_default, "negative", "auto", "no --gtf")]
    sp = spans()

    def body_ok(c, s, e):
        lo, hi = s - NEG_FLANK_FREE, e + NEG_FLANK_FREE
        return peak_mask.covered(c, lo, hi) == 0 and blacklist.covered(c, lo, hi) == 0

    def gene_window(c, s, e):
        return centred(c, s, e, sizes, width=max(NEG_MIN_WIDTH, win_width(e - s)))

    sig = signal() if bw_rows else None
    if all(is_repressive(m) for m in mark_order):
        def hk_gen():
            for gname in HOUSEKEEPING:
                try:
                    g = lookup_gene(sp, gname)
                except LocusError:
                    g = None
                if g is None or not allowed(g[0]) or not body_ok(*g):
                    continue
                w = gene_window(*g)
                if sig is not None:
                    # folds on the gene body +- NEG_FLANK_FREE (the peak-free extent); whole window kept for reference
                    sm = sig.summary(sig.window(g[0], g[1] - NEG_FLANK_FREE, g[2] + NEG_FLANK_FREE))
                    if not sm["mark_fold"] or max(sm["mark_fold"].values()) > NEG_HK_MAX_FOLD:
                        continue
                    kv = sig.kv(sm)
                    kv["worst_mark_fold"] = max(sm["mark_fold"].values())
                    kv["win_target_fold"] = sig.summary(sig.window(*w))["target_fold"]
                else:
                    kv = {}
                yield w, {"name": "negative_" + gname, "anchor": "%s:%d-%d" % g, "anchor_name": gname,
                          "kv": dict(kv, mode="housekeeping", gene_len=g[2] - g[1])}
        rec = place("negative_housekeeping", "negative", "auto", hk_gen(), chooser)
        if rec["window"] is not None:
            return [rec]
        log("negative: no housekeeping gene qualifies (%s); data-driven" % rec["status"])

    cands = []
    for gname, (chs, s, e) in sp.items():
        if len(chs) != 1 or e - s < NEG_MIN_GENE:
            continue
        c = next(iter(chs))
        if not allowed(c) or not body_ok(c, s, e):
            continue
        cands.append((gname, c, s, e))
    if sig is not None:
        scored = []
        for gname, c, s, e in cands:
            # gates and folds on the gene body +- NEG_FLANK_FREE (the peak-free extent); whole window kept for reference
            sm = sig.summary(sig.window(c, s - NEG_FLANK_FREE, e + NEG_FLANK_FREE))
            if sm["igg_fold"] is not None:
                if sm["igg_fold"] < NEG_IGG_MIN_FOLD:
                    continue
            elif sm["rel_cov_max"] is None or sm["rel_cov_max"] < NEG_MIN_REL_COVER:
                continue  # no IgG: some target library must read the region at >= half its typical breadth
            if not sm["mark_fold"]:
                continue
            worst = max(sm["mark_fold"].values())
            kv = sig.kv(sm)
            kv["worst_mark_fold"] = worst
            kv["win_target_fold"] = sig.summary(sig.window(*gene_window(c, s, e)))["target_fold"]
            scored.append((worst, gname, s, c, e, kv))
        scored.sort(key=lambda t: (t[0], t[1], t[2]))
        ordered = [(t[1], t[3], t[2], t[4], t[5]) for t in scored]
        mode = "data-driven lowest worst-mark fold"
    else:
        ordered = sorted(cands, key=lambda t: (-(t[3] - t[2]), t[0], t[2]))
        ordered = [(t[0], t[1], t[2], t[3], {}) for t in ordered]
        mode = "data-driven longest gene"

    def dd_gen():
        for gname, c, s, e, kv in ordered:
            yield gene_window(c, s, e), {"name": "negative_" + gname, "anchor": "%s:%d-%d" % (c, s, e),
                                         "anchor_name": gname, "kv": dict(kv, mode=mode, gene_len=e - s)}

    rec = place("negative_" + "data_driven", "negative", "auto", dd_gen(), chooser)
    if rec["window"] is None:
        log("negative: no data-driven gene passed the filters (%s)" % rec["status"])
    return [rec]


def write_outputs(args, recs, params, log):
    """Cap (SPEC 4.8), numbering and the five output files."""
    os.makedirs(args.outdir, exist_ok=True)
    n_user = sum(1 for r in recs if r["window"] is not None and r.get("user"))
    room = max(0, args.max_regions - n_user)
    final = []
    for r in recs:
        if r["window"] is None:
            continue
        if r.get("user"):
            final.append(r)
        elif room > 0:
            room -= 1
            final.append(r)
        else:
            r["status"] = "skipped:max_regions=%d" % args.max_regions
            r["window"] = None
            log.append("cap: %s dropped (--max-regions %d)" % (r["name"], args.max_regions))
    kept = set(id(r) for r in final)
    for r in recs:
        if r["window"] is None:
            log.append("%s: %s" % (r["name"], r["status"]))
    with open(os.path.join(args.outdir, "regions.tsv"), "w") as fh:
        for k, r in enumerate(final, 1):
            fh.write("%02d_%s\t%s\n" % (k, r["name"], fmt_window(r["window"])))
    with open(os.path.join(args.outdir, "tracks.loci.txt"), "w") as fh:
        fh.write(",".join(fmt_window(r["window"]) for r in final) + "\n")
    with open(os.path.join(args.outdir, "regions.detail.tsv"), "w") as fh:
        fh.write("name\tcategory\tsource\tanchor\tanchor_name\tkey_values\trank_tried\tn_rejected\tstatus\n")
        for r in recs:
            status = "ok" if id(r) in kept else r["status"]
            kv = ";".join("%s=%s" % (k, fmt(v)) for k, v in r["kv"].items())
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n" % (
                r["name"], r["category"], r["source"], r["anchor"], r["anchor_name"], kv, r["rank"], r["rej"], status))
    params["n_regions"] = len(final)
    with open(os.path.join(args.outdir, "params.json"), "w") as fh:
        json.dump(params, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(os.path.join(args.outdir, "peak_qc_regions.log"), "w") as fh:
        fh.write("\n".join(log) + "\n")


def main(argv=None):
    args = parse_args(argv)
    log = []
    try:
        os.makedirs(args.outdir, exist_ok=True)
        inputs = [p for p in (args.peaksets, args.classes, args.bigwigs, args.gtf, args.chrom_sizes, args.config) if p]
        recs, loci_source, loci, neg, neg_source, ref, shape = run(args, log.append, inputs)
        md5 = {}
        for path in inputs:
            md5[path] = md5_of(path)
        params = {
            "argv": sys.argv[1:] if argv is None else list(argv),
            "peaksets": args.peaksets, "classes": args.classes, "bigwigs": args.bigwigs, "gtf": args.gtf,
            "chrom_sizes": args.chrom_sizes, "blacklist": split_list(args.blacklist),
            "exclude_bed": split_list(args.exclude_bed), "nopeak_igg_max_fold": args.nopeak_igg_max_fold,
            "exclude_chrom_regex": args.exclude_chrom_regex, "spike_prefix": args.spike_prefix or None,
            "mito": args.mito or None, "class_names": split_list(args.class_names),
            "shape_marks": shape, "reference_sets": dict((m, ref[m]["id"]) for m in ref),
            "loci_source": loci_source,
            "loci": [{"name": l["name"], "spec": l["spec"], "coord": l["coord"], "flank": l["flank"]} for l in loci],
            "negative": neg if isinstance(neg, str) else {"name": neg["name"], "spec": neg["spec"],
                                                          "flank": neg["flank"]},
            "negative_source": neg_source, "max_regions": args.max_regions, "outdir": args.outdir,
            "constants": {"ROUND_BP": ROUND_BP, "WIDTH_TIERS": WIDTH_TIERS, "WIDTH_STEP": WIDTH_STEP,
                          "DENSEST_WIDTH": DENSEST_WIDTH, "OVERLAP_MAX": OVERLAP_MAX,
                          "BLACKLIST_MAX": BLACKLIST_MAX, "DOMINANT_FRAC": DOMINANT_FRAC,
                          "CROSS_MIN_OVERLAP": CROSS_MIN_OVERLAP, "CROSS_MIN_PAIRS": CROSS_MIN_PAIRS,
                          "NOPEAK_TILE": NOPEAK_TILE, "NOPEAK_STEP": NOPEAK_STEP, "NOPEAK_PAD": NOPEAK_PAD,
                          "NOPEAK_MIN_REL_COVER": NOPEAK_MIN_REL_COVER, "MAX_TRIES": MAX_TRIES,
                          "MIN_PEAKS_SHAPE": MIN_PEAKS_SHAPE, "NEG_MIN_GENE": NEG_MIN_GENE,
                          "NEG_FLANK_FREE": NEG_FLANK_FREE, "NEG_MIN_WIDTH": NEG_MIN_WIDTH,
                          "NEG_IGG_MIN_FOLD": NEG_IGG_MIN_FOLD, "NEG_MIN_REL_COVER": NEG_MIN_REL_COVER,
                          "NEG_HK_MAX_FOLD": NEG_HK_MAX_FOLD, "LOCUS_FLANK": LOCUS_FLANK,
                          "HOUSEKEEPING": HOUSEKEEPING, "POLYCOMB_PATTERNS": POLYCOMB_PATTERNS,
                          "REPRESSIVE_PATTERNS": REPRESSIVE_PATTERNS, "DEFAULT_LOCI": DEFAULT_LOCI},
            "md5": md5,
        }
        write_outputs(args, recs, params, log)
    except PeakQCExit as ex:
        return ex.code
    return 0


if __name__ == "__main__":
    sys.exit(main())
