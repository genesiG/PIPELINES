#!/usr/bin/env python3
r"""make_tracks_ini.py -- pyGenomeTracks tracks.ini and locus regions for one track set.

Called by process TRACKS_INI (tracks.csv row = one call):
    make_tracks_ini.py --track-set NAME --loci "GENEA,GENEB,chr1:1000-2000" --flanks "125000,25000,0" \
        --gtf genes.gtf --gtf-name genes.gtf --samples id1,id2 --bigwigs f1.bw,f2.bw \
        --labels "a,b" --colors "#hex,#hex" [--max-values "100,group:shared_K27me3"] [--bed-files tracks_bed_1/x.bed] \
        [--bed-colors "#hex"] [--bed-labels "title"] [--scalebar BP|auto] [--log2-bigwigs name.bw --log2-for id1] \
        [--number-of-bins 700] [--height 1.0] [--fontsize 14] [--bed-height 0.1] [--scalebar-fontsize 16] \
        [--genes-height 3] [--genes-fontsize 15] [--log2-min -6] [--log2-max 2] [--genes-style UCSC] \
        [--merge-transcripts true] [--gene-rows N] [--arrow-interval N] [--arrowhead-included false] \
        [--bed-after SAMPLE_ID] [--gap-after "ID=H,ID=H"]

Arguments:
  --track-set         track set NAME; prefix of both outputs
  --loci              comma list of exact gene_name values and/or regions chr:start-end (commas allowed inside numbers);
                      an entry may be prefixed LABEL= (e.g. 03_K27me3_largest=chr1:1000-5000): the label becomes the
                      locus slug, so the file names carry it; labels match [A-Za-z0-9_.-]+ and are unique
  --flanks            bp added on each side of a gene; one value or one per locus; ignored for chr:start-end
  --gtf               GTF read for gene coordinates (any path)
  --gtf-name          the same GTF as staged next to the ini (relative); the file of the [Genes] section
  --samples           sample ids, comma list
  --bigwigs           STAGED bigWig file names (relative), one per sample, same order as --samples
  --labels            display label per sample (title of its bigWig track)
  --colors            hex colour per sample
  --max-values        optional; empty = no max_value line; one entry or one per sample. Each entry is one of:
                      number (fixed max_value); auto or empty (no max_value line = pyGenomeTracks autoscale);
                      group:NAME (written as 'max_value = group:NAME'; tracks_limits.py later replaces it by the
                      shared maximum of the group at the plotted locus). NAME matches [A-Za-z0-9_.-]+.
  --bed-files         optional comma list of STAGED BED files (relative); one bed track each
  --bed-colors        optional comma list, one hex colour per --bed-files entry (used for color and border_color);
                      default #8c96c6 for every bed
  --bed-labels        optional comma list, one title per --bed-files entry; default = file stem
  --scalebar          optional scale bar size in bp, or auto (pyGenomeTracks picks the size: no size line);
                      empty = no [scalebar] section
  --log2-bigwigs      optional comma list of STAGED log2 bigWig names
  --log2-for          optional comma list of sample ids, one per --log2-bigwigs entry: that log2 track is drawn after this sample
  --number-of-bins    positive integer N, written as number_of_bins = N in every sample bigWig section (default 700)
  --height            height of each sample bigWig track (default 1.0)
  --fontsize          accepted for compatibility; NOT written (pyGenomeTracks bigWig tracks ignore it; their title
                      size is the global pyGenomeTracks --fontSize) (default 14)
  --bed-height        height of each bed track (default 0.1)
  --scalebar-fontsize font size of the scale bar (default 16)
  --genes-height      height of the [Genes] track (default 3)
  --genes-fontsize    font size of the [Genes] track labels (default 15)
  --genes-style       [Genes] style: UCSC (default), flybase or tssarrow
  --merge-transcripts true (default, one model per gene) or false (every transcript is a model)
  --gene-rows         optional; written as gene_rows = N (caps the rows of the stacked gene display; genes that do not fit are hidden)
  --arrowhead-included  true or false (default false); true writes arrowhead_included = true
  --arrow-interval    optional; written as arrow_interval = N in [Genes] (pyGenomeTracks default 2; larger = fewer intron arrows)
  --bed-after         optional sample id: the BED track(s) are drawn right after this sample's track (0.15 spacer before and
                      after each BED) instead of after the last sample
  --gap-after         optional "ID=HEIGHT,ID=HEIGHT" (comma or semicolon): spacer height after that sample's track (default 0.25); one
                      entry per antibody boundary reproduces the legacy 0.5 group divider
  --log2-min          min_value of the log2 tracks (default -6); log2 tracks keep their own fixed limits, never grouped
  --log2-max          max_value of the log2 tracks (default 2)

Gene lookup: a GTF line belongs to gene X iff its attribute column (9th field) contains exactly the token
  gene_name "X";   (regex (?:^|[;\s])gene_name "X";, X escaped). The region is the min start .. max end over those
  lines (any feature type). The chromosome must be unique. Region = max(1, start - flank) .. end + flank.
  Unknown gene, or a gene on several chromosomes -> exit 3.
  Locus chr:start-end is used as given.

Outputs (working directory):
  NAME.tracks.ini   relative file names only (no absolute paths). Section order: [scalebar] (if given; for auto no size
                    line) and a spacer;
                    per sample: bigwig section (title = label, colour, min_value 0, max_value per the
                    --max-values entry (no line for auto or empty), number_of_bins, height, nans_to_zeros,
                    show_data_range true; no fontsize line, see --fontsize); if the sample has a log2 bigWig, a log2
                    section (color red, negative_color blue, min_value --log2-min, max_value --log2-max, show_data_range
                    true) followed by a zero-line section (file_type hlines, overlay_previous share-y, show_data_range
                    false); then a spacer. Then one bed section per BED (type bed, display collapsed,
                    labels no, style UCSC, height --bed-height, color and border_color from --bed-colors, titled from
                    --bed-labels or by the file stem) with a spacer after each, a spacer, and [Genes] (type genes, file =
                    --gtf-name, labels yes, prefered_name gene_name, display stacked, merge_transcripts true (one model per
                    gene; no gene_rows cap, which hides genes), style UCSC, fontstyle italic, fontsize --genes-fontsize,
                    height --genes-height).
                    No [x-axis]. All section names are unique.
  NAME.regions.tsv  one line per locus: NAME.<locus_slug><TAB>chr:start-end.
                    locus_slug = <gene>_<flank/1000>kb if flank % 1000 == 0, else <gene>_<flank>bp;
                    for chr:start-end loci: chr_start_end (commas removed, ':' and '-' -> '_');
                    an entry prefixed LABEL= uses LABEL as its slug instead.
Stdout: one line per locus, "<locus> -> <chr:start-end>".
Exit 3: unknown gene; gene on several chromosomes; a list length that does not match (labels, colors, bigwigs = number
        of samples; max-values and flanks = 1 or the required count; bed-colors and bed-labels = number of bed files);
        an empty bed colour or label entry; a max-values entry that is not a number, auto, empty or group:NAME; a group
        name outside [A-Za-z0-9_.-]+; a numeric option (height, fontsize, bed-height, scalebar-fontsize, genes-height,
        genes-fontsize, log2-min, log2-max) that is not a number; number-of-bins that is not a positive integer;
        non-integer flank;
        malformed locus; malformed or duplicate locus label; absolute file name; GTF missing; duplicate sample id;
        log2 entries that do not match the samples.
Exit 2: a required argument is missing.
"""
import argparse
import os
import re
import sys

DEFAULT_BED_COLOR = "#8c96c6"
LOG2_COLOR = "#e34948"
LOG2_NEG_COLOR = "#2a78d6"
ZERO_LINE_COLOR = "#8a8984"
GENE_COLOR = "#000000"
ATTR_TOKEN = r'(?:^|[;\s])gene_name "%s";'
LOCUS_RE = re.compile(r"^([^:\s]+):([\d,]+)-([\d,]+)$")
GROUP_NAME_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def split_plain(s):
    if s is None or s.strip() == "":
        return []
    return [t.strip() for t in s.split(",")]


def split_loci(s):
    # A comma that is a thousands separator (followed by three digits and then ',', '-' or the end) is kept.
    if s is None or s.strip() == "":
        return []
    return [t.strip() for t in re.split(r",(?!\d{3}(?:,|-|$))", s)]


def to_int(text, tok):
    try:
        return int(text.replace(",", ""))
    except ValueError:
        die("locus '%s' has a non-integer coordinate" % tok)


def expand(opt, values, n):
    # one value is used for every locus or sample; otherwise one value per entry is required
    if len(values) == 1:
        return values * n
    if len(values) != n:
        die("%s has %d entries; expected 1 or %d" % (opt, len(values), n))
    return list(values)


def need_per_sample(opt, values, n):
    if len(values) != n:
        die("%s has %d entries but --samples has %d" % (opt, len(values), n))
    if any(v == "" for v in values):
        die("%s has an empty entry" % opt)


def need_per_bed(opt, values, n):
    # no values = option omitted (the caller supplies its default)
    if not values:
        return []
    if len(values) != n:
        die("%s has %d entries but --bed-files has %d" % (opt, len(values), n))
    if any(v == "" for v in values):
        die("%s has an empty entry" % opt)
    return list(values)


def check_relative(opt, path):
    if path == "" or os.path.isabs(path):
        die("%s must be a relative file name (got '%s')" % (opt, path))


def check_number(opt, text):
    try:
        float(text)
    except ValueError:
        die("%s is not a number (got '%s')" % (opt, text))


def max_value_token(v):
    """One --max-values entry: the text of a number or group:NAME, or None (no max_value line = autoscale)."""
    if v in ("", "auto", "empty"):
        return None
    if v.startswith("group:"):
        if not GROUP_NAME_RE.match(v[len("group:"):]):
            die("--max-values group name must match [A-Za-z0-9_.-]+ (got '%s')" % v)
        return v
    try:
        float(v)
    except ValueError:
        die("--max-values entry is not a number, auto, empty or group:NAME (got '%s')" % v)
    return v


def gene_regions(gtf, genes):
    """Per gene: set of chromosomes and min start / max end over the GTF lines carrying the exact gene_name token."""
    patterns = {g: re.compile(ATTR_TOKEN % re.escape(g)) for g in genes}
    needles = {g: 'gene_name "%s";' % g for g in genes}
    found = {g: {"chroms": set(), "start": None, "end": None} for g in genes}
    with open(gtf, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 9:
                continue
            attr = f[8]
            if 'gene_name "' not in attr:
                continue
            for g in genes:
                if needles[g] not in attr or not patterns[g].search(attr):
                    continue
                try:
                    s, e = int(f[3]), int(f[4])
                except ValueError:
                    die("non-integer start/end on GTF line of gene %s: %s" % (g, line.rstrip()))
                rec = found[g]
                rec["chroms"].add(f[0])
                rec["start"] = s if rec["start"] is None else min(rec["start"], s)
                rec["end"] = e if rec["end"] is None else max(rec["end"], e)
    return found


def render(name, items, names):
    names.append(name)
    lines = ["[%s]" % name]
    for k, v in items:
        lines.append("%s = %s" % (k, v) if v != "" else "%s =" % k)
    return "\n".join(lines)


def build_ini(samples, labels, bigwigs, colors, max_vals, log2_of, beds, bed_colors, bed_labels, scalebar, gtf_name, o,
              bed_after="", gap_after=None):
    gap_after = gap_after or {}
    blocks, names, nsp = [], [], [0]

    def spacer(h):
        nsp[0] += 1
        blocks.append(render("spacer %d" % nsp[0], [("file_type", "spacer"), ("height", h)], names))

    if scalebar:
        sb_items = [("file_type", "scalebar"), ("where", "top")]
        if scalebar != "auto":
            sb_items.append(("size", scalebar))
        sb_items += [
            ("line_width", 2), ("color", "black"), ("fontsize", o.scalebar_fontsize),
            ("height", 0.9), ("title", ""),
        ]
        blocks.append(render("scalebar", sb_items, names))
        spacer(0.25)

    def emit_beds():
        for k, b in enumerate(beds):
            stem = os.path.splitext(os.path.basename(b))[0]
            # pyGenomeTracks guesses the file type from the extension and only knows .bed; other BED-like files
            # (MACS3 .broadPeak / .narrowPeak, ...) need the type given explicitly
            file_type = [] if os.path.splitext(b)[1].lower() == ".bed" else [("file_type", "bed")]
            blocks.append(render("bed %d %s" % (k + 1, stem), [
                ("file", b)] + file_type + [
                ("title", bed_labels[k]), ("fontsize", 14), ("color", bed_colors[k]),
                ("border_color", bed_colors[k]), ("type", "bed"), ("height", o.bed_height),
                ("display", "collapsed"), ("labels", "no"), ("style", "UCSC"),
            ], names))
            spacer(0.15)

    for i, sid in enumerate(samples):
        items = [("file", bigwigs[i]), ("title", labels[i]), ("color", colors[i]),
                 ("type", "fill"), ("height", o.height), ("number_of_bins", o.number_of_bins), ("min_value", 0)]
        if max_vals and max_vals[i] is not None:
            items.append(("max_value", max_vals[i]))
        items += [("show_data_range", "true"), ("nans_to_zeros", "true")]
        blocks.append(render(sid, items, names))
        if sid in log2_of:
            blocks.append(render("%s log2" % sid, [
                ("file", log2_of[sid]), ("title", "%s log2 ratio" % labels[i]),
                ("color", LOG2_COLOR), ("negative_color", LOG2_NEG_COLOR), ("type", "fill"),
                ("height", 1.0), ("min_value", o.log2_min), ("max_value", o.log2_max),
                ("show_data_range", "true"), ("nans_to_zeros", "true"),
            ], names))
            blocks.append(render("%s log2 zero line" % sid, [
                ("file_type", "hlines"), ("y_values", 0), ("overlay_previous", "share-y"),
                ("color", ZERO_LINE_COLOR), ("line_width", 0.75), ("show_data_range", "false"),
            ], names))
        if bed_after and sid == bed_after:
            spacer(0.15)
            emit_beds()
        else:
            spacer(gap_after.get(sid, 0.25))

    if not bed_after:
        emit_beds()

    spacer(0.5)
    gene_items = [
        ("file", gtf_name), ("title", ""), ("color", GENE_COLOR), ("height", o.genes_height), ("type", "genes"),
        ("labels", "yes"), ("fontsize", o.genes_fontsize), ("prefered_name", "gene_name"),
    ]
    if o.gene_rows:
        gene_items.append(("gene_rows", o.gene_rows))
    gene_items += [("display", "stacked"), ("merge_transcripts", o.merge_transcripts),
                   ("style", o.genes_style), ("fontstyle", "italic")]
    if o.arrow_interval:
        gene_items.append(("arrow_interval", o.arrow_interval))
    if o.arrowhead_included == "true":
        gene_items.append(("arrowhead_included", "true"))
    blocks.append(render("Genes", gene_items, names))

    if len(set(names)) != len(names):
        die("internal: duplicate section names in the ini")
    return "\n\n".join(blocks) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--track-set", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--flanks", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--gtf-name", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--bigwigs", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--colors", required=True)
    ap.add_argument("--max-values", default="")
    ap.add_argument("--bed-files", default="")
    ap.add_argument("--bed-colors", default="")
    ap.add_argument("--bed-labels", default="")
    ap.add_argument("--scalebar", default="")
    ap.add_argument("--log2-bigwigs", default="")
    ap.add_argument("--log2-for", default="")
    ap.add_argument("--number-of-bins", default="700")
    ap.add_argument("--height", default="1.0")
    ap.add_argument("--fontsize", default="14")
    ap.add_argument("--bed-height", default="0.1")
    ap.add_argument("--scalebar-fontsize", default="16")
    ap.add_argument("--genes-height", default="3")
    ap.add_argument("--genes-fontsize", default="15")
    ap.add_argument("--genes-style", default="UCSC", choices=["UCSC", "flybase", "tssarrow"])
    ap.add_argument("--merge-transcripts", default="true", choices=["true", "false"])
    ap.add_argument("--gene-rows", default="")
    ap.add_argument("--arrowhead-included", default="false", choices=["true", "false"])
    ap.add_argument("--arrow-interval", default="")
    ap.add_argument("--bed-after", default="")
    ap.add_argument("--gap-after", default="")
    ap.add_argument("--log2-min", default="-6")
    ap.add_argument("--log2-max", default="2")
    a = ap.parse_args()

    if not os.path.isfile(a.gtf):
        die("GTF not found: %s" % a.gtf)

    # samples, labels, colours, bigWigs
    samples = split_plain(a.samples)
    n = len(samples)
    if n == 0 or any(s == "" for s in samples):
        die("--samples is empty or has an empty id")
    if len(set(samples)) != n:
        die("--samples has duplicate ids")
    labels = split_plain(a.labels)
    colors = split_plain(a.colors)
    bigwigs = split_plain(a.bigwigs)
    need_per_sample("--labels", labels, n)
    need_per_sample("--colors", colors, n)
    need_per_sample("--bigwigs", bigwigs, n)
    for b in bigwigs:
        check_relative("--bigwigs", b)
    max_vals = split_plain(a.max_values)
    if max_vals:
        max_vals = [max_value_token(v) for v in expand("--max-values", max_vals, n)]

    # styling options
    for opt, val in (("--height", a.height), ("--fontsize", a.fontsize), ("--bed-height", a.bed_height),
                     ("--scalebar-fontsize", a.scalebar_fontsize), ("--genes-height", a.genes_height),
                     ("--genes-fontsize", a.genes_fontsize), ("--log2-min", a.log2_min), ("--log2-max", a.log2_max)):
        check_number(opt, val)
    if not a.number_of_bins.isdigit() or int(a.number_of_bins) <= 0:
        die("--number-of-bins must be a positive integer (got '%s')" % a.number_of_bins)

    # beds, log2 tracks, gtf name, scale bar
    beds = split_plain(a.bed_files)
    for b in beds:
        check_relative("--bed-files", b)
    bed_colors = need_per_bed("--bed-colors", split_plain(a.bed_colors), len(beds)) or [DEFAULT_BED_COLOR] * len(beds)
    bed_labels = need_per_bed("--bed-labels", split_plain(a.bed_labels), len(beds)) or \
        [os.path.splitext(os.path.basename(b))[0] for b in beds]
    check_relative("--gtf-name", a.gtf_name)
    log2_bws = split_plain(a.log2_bigwigs)
    log2_for = split_plain(a.log2_for)
    if len(log2_bws) != len(log2_for):
        die("--log2-bigwigs has %d entries but --log2-for has %d" % (len(log2_bws), len(log2_for)))
    for b in log2_bws:
        check_relative("--log2-bigwigs", b)
    for sid in log2_for:
        if sid not in samples:
            die("--log2-for '%s' is not one of --samples" % sid)
    if len(set(log2_for)) != len(log2_for):
        die("--log2-for lists a sample more than once")
    log2_of = dict(zip(log2_for, log2_bws))
    scalebar = a.scalebar.strip()
    if scalebar and scalebar != "auto" and (not scalebar.isdigit() or int(scalebar) <= 0):
        die("--scalebar must be auto or a positive bp value (got '%s')" % scalebar)

    # loci and flanks
    loci = split_loci(a.loci)
    if not loci or any(t == "" for t in loci):
        die("--loci is empty or has an empty entry")
    # optional LABEL= prefix: the label is the locus slug; the rest of the code sees the stripped token
    loci_labels = []
    for i, tok in enumerate(loci):
        label = ""
        if "=" in tok:
            label, tok = (s.strip() for s in tok.split("=", 1))
            if not GROUP_NAME_RE.match(label):
                die("locus label '%s' must match [A-Za-z0-9_.-]+" % label)
            if tok == "":
                die("--loci is empty or has an empty entry")
            if label in loci_labels:
                die("duplicate locus label '%s'" % label)
            loci[i] = tok
        loci_labels.append(label)
    flanks = split_plain(a.flanks)
    for f in flanks:
        if not f.isdigit():
            die("flank must be a non-negative integer bp (got '%s')" % f)
    flanks = expand("--flanks (one per locus)", flanks, len(loci))

    parsed = {}
    gene_set = []
    for tok in loci:
        if ":" in tok:
            m = LOCUS_RE.match(tok)
            if not m:
                die("locus '%s' is neither a gene name nor chr:start-end" % tok)
            start, end = to_int(m.group(2), tok), to_int(m.group(3), tok)
            if start < 1 or end <= start:
                die("locus '%s': need 1 <= start < end" % tok)
            parsed[tok] = (m.group(1), start, end)
        elif tok not in gene_set:
            gene_set.append(tok)

    found = gene_regions(a.gtf, gene_set) if gene_set else {}
    unknown = [g for g in gene_set if not found[g]["chroms"]]
    if unknown:
        die("gene(s) not found in %s (exact gene_name): %s" % (a.gtf, ", ".join(unknown)))
    for g in gene_set:
        if len(found[g]["chroms"]) != 1:
            die("gene %s lies on several chromosomes (%s); cannot build one region"
                % (g, ", ".join(sorted(found[g]["chroms"]))))

    records = []  # (locus token, locus slug, region)
    for tok, fl, label in zip(loci, flanks, loci_labels):
        if ":" in tok:
            chrom, start, end = parsed[tok]
            slug = tok.replace(",", "").replace(":", "_").replace("-", "_")
        else:
            rec = found[tok]
            chrom = next(iter(rec["chroms"]))
            fl_int = int(fl)
            start = max(1, rec["start"] - fl_int)
            end = rec["end"] + fl_int
            if fl_int % 1000 == 0:
                slug = "%s_%dkb" % (tok, fl_int // 1000)
            else:
                slug = "%s_%dbp" % (tok, fl_int)
        if label:
            slug = label
        records.append((tok, slug, "%s:%d-%d" % (chrom, start, end)))

    if a.gene_rows and not (a.gene_rows.isdigit() and int(a.gene_rows) > 0):
        die("--gene-rows must be a positive integer (got '%s')" % a.gene_rows)
    if a.bed_after and a.bed_after not in samples:
        die("--bed-after '%s' is not one of --samples" % a.bed_after)
    if a.bed_after and not beds:
        die("--bed-after needs --bed-files")
    gap_after = {}
    for ent in [t.strip() for t in re.split(r"[;,]", a.gap_after) if t.strip()]:
        sid, _, h = ent.partition("=")
        if sid.strip() not in samples:
            die("--gap-after '%s' is not one of --samples" % sid.strip())
        check_number("--gap-after", h.strip())
        gap_after[sid.strip()] = h.strip()
    text = build_ini(samples, labels, bigwigs, colors, max_vals, log2_of, beds, bed_colors, bed_labels,
                     scalebar, a.gtf_name, a, a.bed_after, gap_after)
    with open(a.track_set + ".tracks.ini", "w", encoding="utf-8") as fh:
        fh.write(text)
    with open(a.track_set + ".regions.tsv", "w", encoding="utf-8") as fh:
        for _, slug, region in records:
            fh.write("%s.%s\t%s\n" % (a.track_set, slug, region))
    for tok, _, region in records:
        print("%s -> %s" % (tok, region))


if __name__ == "__main__":
    main()
