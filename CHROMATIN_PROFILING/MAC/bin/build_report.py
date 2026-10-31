#!/usr/bin/env python3
"""build_report.py -- assemble the self-contained HTML QC report from the QC tables and SVG plots.

Called by REPORT (modules/local/report.nf):
    build_report.py --project NAME --tables-dir . --plots-dir . --outdir .

Arguments:
  --project     project name (page title and output file name)
  --tables-dir  directory holding the QC tables (*.tsv): alignment_stats, fragment_sizes, library_complexity, frip,
                peak_stats, tss_enrichment, correlation_matrix. Missing tables are omitted, not an error.
  --plots-dir   directory holding the QC plots (*.svg). Missing plots are omitted.
  --outdir      where to write the report (default ".")
Output: <project>_qc_report.html, one self-contained file (python3 standard library only).
  - every table is a plain HTML <table>;
  - every SVG is inlined as raw <svg>. Each plot gets a unique id prefix: every id="X" becomes id="<prefix>-X" and
    every href="#X", url(#X) and CSS #X selector that points to it is rewritten (several SVGs in one page would
    otherwise reuse short ids such as glyph0-0 or clip1 and the browser would mix their shapes);
  - every kind of plot has a text block: what it measures, what the ideal case looks like, what a deviation means
    (bullet list), and edge cases. Plots of one kind (per-group fingerprints, per-sample TSS profiles) share one block.
Exit 3 when neither a table nor a plot is found.
"""
import argparse
import html
import os
import re
import sys

# --------------------------------------------------------------------------------------------------------------
# Text blocks. Style (porting guide section 7): direct declarative sentences; what it measures, then the ideal case,
# then what a deviation means; troubleshooting as bullets; edge cases covered; no comparison language unless the
# plot is about a comparison.
# --------------------------------------------------------------------------------------------------------------
PLOTS = {
    "alignment_read_composition": dict(
        title="Read composition (counts)",
        measures="Each bar splits the raw reads of one sample, in millions, into four classes: target-genome reads, "
                 "mitochondrial reads (chrM), spike-in reads, and unmapped reads. Raw means counted on the aligned BAM "
                 "before MAPQ filtering, duplicate removal and chrM/spike-in removal.",
        ideal="Most of each bar is target-genome reads. Unmapped and chrM reads are small slices. The spike-in slice is "
              "present and sized as designed.",
        deviations=["A large unmapped slice means many reads did not align to the combined index: check adapter trimming, "
                    "the index and the species of the sample.",
                    "A large chrM slice means mitochondrial DNA took up sequencing depth; the usable target reads are "
                    "fewer than the bar height suggests.",
                    "A spike-in slice that is much larger or smaller in one sample than in the others of the same antibody "
                    "points to a spike-in addition problem (wrong amount, missed addition) rather than biology.",
                    "A short bar means a shallow library; read counts feed every later step, so very short bars limit "
                    "peak calling and normalisation."],
        edge=["Samples built from pre-aligned BAMs have no raw statistics and do not appear in this plot.",
              "Without a spike-in contig prefix the spike-in slice is empty for every sample."]),
    "alignment_read_composition_pct": dict(
        title="Read composition (percent)",
        measures="The same four classes as the count plot, each bar scaled to 100% of the raw reads of the sample.",
        ideal="A dominant target-genome share with small unmapped and chrM shares, and a spike-in share that is similar "
              "between samples of the same antibody.",
        deviations=["The target share falls when unmapped, chrM or spike-in reads grow; read the cause from the "
                    "slice that changed.",
                    "A spike-in share that varies strongly within one antibody reflects different spike-in to sample "
                    "ratios and affects spike-in normalisation."],
        edge=["A sample with very few reads still fills 100%; read this plot together with the count plot."]),
    "alignment_pct_spikein": dict(
        title="Spike-in reads per sample, by antibody",
        measures="The percentage of mapped reads that fall on spike-in contigs, one panel per antibody. The spike-in "
                 "percentage depends on how much carrier DNA was added relative to the target chromatin, so each "
                 "antibody has its own typical level.",
        ideal="Samples within a panel (same antibody) show similar percentages. Replicates of one condition sit close "
              "together.",
        deviations=["One sample far above its panel mates has less target signal or more spike-in than intended.",
                    "One sample far below its panel mates has more target signal or less spike-in than intended.",
                    "A control antibody (IgG) usually has a higher spike-in percentage than a specific antibody because "
                    "it captures little target chromatin; that difference between panels is expected."],
        edge=["A single sample in a panel has no within-antibody reference; judge it against the read-count plot.",
              "Percentages near zero mean spike-in counts are too low for spike-in normalisation (see the count plot)."]),
    "alignment_n_spikein": dict(
        title="Spike-in read counts",
        measures="The number of mapped reads on spike-in contigs per sample. The dashed line marks 1,000 reads.",
        ideal="Every sample has well over 1,000 spike-in reads, so spike-in based scale factors rest on enough counts.",
        deviations=["Below 1,000 reads the spike-in count is noisy. Scale factors derived from it are unreliable.",
                    "Zero reads means no spike-in was added, the spike-in contigs are absent from the index, or the contig "
                    "prefix given to the pipeline does not match the index."],
        edge=["Counts are raw (before deduplication); a library with heavy duplication can show a high count from few "
              "distinct molecules."]),
    "alignment_pct_chrM": dict(
        title="Mitochondrial reads",
        measures="The percentage of mapped reads on the mitochondrial contig, per sample.",
        ideal="A low percentage. Nuclear chromatin assays recover little mitochondrial DNA.",
        deviations=["A high percentage in one sample means the nuclear fraction was diluted by mitochondrial DNA; "
                    "usable depth is lower than total depth.",
                    "A high percentage across all samples points to the protocol (cell lysis, nuclear prep) or to the "
                    "assay type."],
        edge=["A value of zero in every sample suggests the mitochondrial contig name in the index differs from the one "
              "given to the pipeline."]),
    "alignment_duplicates": dict(
        title="Duplicate fractions (Picard MarkDuplicates)",
        measures="For each sample, the share of read pairs examined by Picard that are unique, PCR duplicates, or "
                 "optical duplicates.",
        ideal="A mostly unique bar with small duplicate shares. Duplicates in low-input CUT&RUN/CUT&Tag are partly "
              "expected because the fragment start and end are fixed by the enzyme cut.",
        deviations=["A large PCR duplicate share means few distinct molecules were amplified many times: low input or "
                    "over-amplification. See the library complexity figure for the same information per fragment.",
                    "A large optical duplicate share points to a flow-cell or clustering issue.",
                    "An antibody-specific shift (IgG high, specific antibody low, or the reverse) reflects how much "
                    "chromatin each antibody captured."],
        edge=["Samples without Picard metrics do not appear."]),
    "alignment_final_reads": dict(
        title="Reads in the final BAM",
        measures="The number of mapped reads left in the final BAM of each sample after MAPQ filtering, duplicate removal "
                 "(where enabled), chrM and spike-in removal, and blacklist filtering.",
        ideal="Depth that supports the planned analysis, similar between replicates.",
        deviations=["Replicates that differ several-fold in final reads give unequal weight in peak calling and "
                    "normalisation.",
                    "A control library (IgG) much deeper than its targets is not harmful; a very shallow control "
                    "gives a noisy background."],
        edge=["Pre-aligned BAM samples appear here although they have no raw statistics."]),
    "fragment_sizes_density": dict(
        title="Fragment size distributions",
        measures="The distribution of fragment lengths (distance between the outer ends of the two mates) of properly "
                 "paired reads in each final BAM. One line per sample, coloured by antibody. The dashed guides mark "
                 "120, 200, 400 and 600 bp.",
        ideal="A clear peak whose position fits the target: roughly 150-200 bp for nucleosome-associated chromatin "
              "marks, with a ladder of smaller peaks at multiples of about 200 bp for di- and tri-nucleosomes. "
              "Transcription-factor targets add a shorter population below 120 bp. Replicates of one antibody trace the "
              "same curve.",
        deviations=["A broad, featureless distribution means the chromatin was not cleanly digested (over- or "
                    "under-digestion).",
                    "A shift of one replicate relative to its mates points to a different digestion time or enzyme amount.",
                    "A dominant peak far above 400 bp means large, poorly cut fragments.",
                    "A peak below 100 bp that dominates a histone-mark sample suggests over-digestion or "
                    "non-nucleosomal background."],
        edge=["IgG controls often show a flatter curve with a broader size range; that is expected for a library "
              "with little specific signal.",
              "Fragments longer than the bowtie2 maximum insert size are not aligned as proper pairs and do not appear "
              "here."]),
    "fragment_sizes_boxplot": dict(
        title="Median fragment size by antibody",
        measures="One point per sample: the median fragment size, grouped by antibody. The box shows the spread of the "
                 "medians within the antibody.",
        ideal="Points of one antibody lie close together, at the size expected for the target.",
        deviations=["A point away from its antibody group marks a sample whose digestion differed; compare with the "
                    "density plot.",
                    "A wide box means the digestion varied across the replicates of that antibody."],
        edge=["With one or two samples per antibody the box has little meaning; read the individual points."]),
    "library_complexity": dict(
        title="Library complexity",
        measures="Four panels, one bar per sample. NRF is the number of distinct fragments divided by all fragments. "
                 "PBC1 is the fraction of distinct fragments seen exactly once. PBC2 is the number of fragments seen "
                 "once divided by the number seen twice. Fragments are defined by chromosome, start and end of "
                 "properly paired read1 alignments before duplicate removal. The last panel shows the duplication "
                 "percentage from Picard EstimateLibraryComplexity. Dashed guides: NRF and PBC1 at 0.8 and 0.9; PBC2 "
                 "at 3 and 10.",
        ideal="NRF and PBC1 above 0.9, PBC2 above 10, and a low Picard duplication percentage. The library holds many "
              "distinct molecules relative to the sequencing depth.",
        deviations=["NRF or PBC1 between 0.8 and 0.9 and PBC2 between 3 and 10 are acceptable for low-input "
                    "libraries.",
                    "Values below those ranges mean most reads come from few molecules: the library is over-amplified "
                    "or the input was very low. More sequencing adds mostly duplicates.",
                    "NRF high but PBC2 low means a few fragment positions are very deep while most are single: look for "
                    "strong pile-ups (for example repeats or a very strong peak)."],
        edge=["These thresholds come from ChIP-seq. Fragment ends in CUT&RUN are enzyme-defined, so strong, narrow "
              "peaks produce identical fragments at high depth and lower NRF without a library problem.",
              "A sample with very few fragments has unstable ratios; PBC2 is undefined (NA) when no fragment occurs "
              "exactly twice.",
              "Panels show 'not available' when Picard output is absent."]),
    "preseq_ccurve": dict(
        title="Preseq complexity curves",
        measures="Distinct fragments as a function of fragments sequenced, one line per sample, from preseq c_curve on "
                 "the filtered BAM before duplicate removal. The dashed diagonal marks a library without any duplicate.",
        ideal="Curves that stay close to the diagonal and keep rising at the end of the data. The library still yields "
              "new molecules with more sequencing.",
        deviations=["A curve that bends flat early means additional sequencing yields few new molecules.",
                    "A curve that leaves the diagonal at the very start means duplication from the first reads on: low "
                    "input or heavy amplification.",
                    "Curves that end at different depths reflect different sequencing depth; compare shapes, not end "
                    "points."],
        edge=["A curve with only a few points belongs to a small library; the step size was reduced automatically "
              "for such libraries.",
              "Strong, narrow peaks make a curve flatter even in a good library because many fragments share ends."]),
    "fingerprint": dict(
        title="Fingerprint (signal concentration)",
        measures="For one group of samples, the genome is cut into 500 bp bins and sorted by read count. The line shows "
                 "the cumulative fraction of reads (y) against the fraction of bins (x), lowest to highest.",
        ideal="A specific antibody bends sharply away from the diagonal: a small fraction of bins holds most of the "
              "reads. The IgG control stays close to the diagonal.",
        deviations=["A specific antibody close to the diagonal means little enrichment: the antibody did not work or "
                    "the sample failed.",
                    "A control that bends strongly means enrichment in the control (open-chromatin bias or "
                    "contamination).",
                    "Broad marks bend gently, narrow marks bend steeply. A gentle bend is the expected shape for a "
                    "broad histone mark."],
        edge=["A group without a control shows the specific samples only; judge enrichment against the diagonal.",
              "A very shallow library draws a jagged line because most bins hold zero or one read."]),
    "correlation_heatmap": dict(
        title="Sample correlation",
        measures="Spearman correlation of read counts in 5 kb bins between every pair of samples, clustered. The numbers "
                 "in the cells are the coefficients. Counts use whole fragments.",
        ideal="Replicates of one antibody and condition cluster together with high correlation. Specific antibodies "
              "separate from the IgG controls.",
        deviations=["A replicate that clusters away from its mates is an outlier: check its other QC metrics.",
                    "Samples cluster by batch or by depth instead of by antibody or condition: a technical effect "
                    "dominates.",
                    "Specific samples correlating strongly with IgG means weak enrichment."],
        edge=["Narrow-peak targets have most bins near zero, so correlations are lower than for broad marks even "
              "between good replicates.",
              "Cells for two nearly empty samples can be high because both are dominated by background."]),
    "tss_profile": dict(
        title="TSS profile",
        measures="Mean normalised signal in a window around annotated transcription start sites (TSS), one plot per "
                 "sample. The enrichment table gives the ratio of the signal at the TSS to the signal in the outer "
                 "flanks.",
        ideal="For a mark found at active promoters, a clear peak at the TSS or a symmetric dip with shoulders, and a "
              "flat low baseline in the flanks. For a mark absent from promoters (for example a repressive mark), a flat "
              "profile is the expected result.",
        deviations=["A flat profile for a promoter-associated mark means no enrichment at TSSs.",
                    "A high flank level compared to the TSS peak means high background.",
                    "An IgG control with a TSS peak means open-chromatin bias."],
        edge=["Scores are only comparable between samples that use the same normalisation and bigWig units.",
              "Repressive and broad domain marks have low scores at TSS by design."]),
    "replicate_correlation": dict(
        title="Replicate concordance",
        measures="For each pair of libraries of the same group (replicates), agreement of read counts in 10 kb bins "
                 "outside the blacklist: Pearson correlation of log2 counts per million, Spearman correlation of the "
                 "counts, the same Pearson on the 10 % of bins with the highest mean signal, and the robust spread of the "
                 "replicate log2 ratio on those bins. Grey points are pairs of libraries of the same antibody from "
                 "different groups.",
        ideal="Replicate pairs agree more closely than pairs from different conditions of the same antibody, and the "
              "scatter is a narrow band along the diagonal.",
        deviations=["A replicate pair no closer than the non-replicate pairs: check depth, fragment sizes and "
                    "fingerprints of both libraries.",
                    "High genome-wide but low signal-bin correlation: agreement comes from where the domains are, not "
                    "from how strong they are."],
        edge=["Values depend on depth, bin size and how broad the mark is; there is no universal cut-off, so compare "
              "within a dataset.",
              "Groups with one library have no value."]),
    "norm_ma": dict(
        title="Normalisation QC (MA plot)",
        measures="For one normalisation run, one MA panel per pair of libraries: each pair of replicates of a group "
                 "(where no biological change is expected) and each treated library against the reference-condition "
                 "library with the same replicate label. Each bin outside the blacklist is a point: x is the mean and y the "
                 "difference (other minus reference, negative = less signal in the other library) of the two libraries' log2 "
                 "counts per million (CPM); the points are drawn as a density. One dashed line per method marks where that "
                 "method declares no change; dotted = factor from a QC-failed or imputed library; numbers in the corner are "
                 "the line positions. Colours: orange for ChIPseqSpikeInFree, blue for TMM, black for CPM, green for spike-in.",
        ideal="After a correct normalisation the background cloud, the dense region at low mean, sits on y = 0, and the "
              "line of the method that was used lies on that cloud. The black bar at the right margin is the measured core "
              "of the cloud, between its density peak and its density-weighted centroid; the white dot marks the density "
              "peak. A method line that lies inside the bar agrees with the cloud.",
        deviations=["A method line outside the cloud means that method puts the two libraries on a different scale from "
                    "most bins. For spike-in and ChIPseqSpikeInFree this can be the correct outcome when the bins in the "
                    "cloud truly changed (for example a global loss of a mark), so read the lines as a comparison of "
                    "methods, not as proof that one is right. For TMM the cloud is centred on zero by construction.",
                    "A cloud centred away from y = 0 for every method points to a difference in background between the two "
                    "libraries that depth does not explain: check the other QC plots of both libraries.",
                    "A dotted line means one of the two libraries failed that method's QC, so its factor is imputed and the "
                    "line is not a measured value.",
                    "A cloud with two separate peaks, or one much wider than the others, mixes bins with and without a real "
                    "change; the core bar then spans a wide range."],
        edge=["Without a pairs table the core bar and the white dot are absent and only the lines are drawn.",
              "Each panel has its own reference library (named in the panel title); panels are not compared with each other.",
              "A library with very few bins above the count filter gives a sparse cloud that says little about the offset."]),
}
GENERIC = dict(
    title="QC plot", measures="A QC plot produced by the pipeline.", ideal="See the table next to it for the values.",
    deviations=["Compare the sample against the other samples of its antibody."], edge=[])

TABLES = {
    "alignment_stats": ("Alignment statistics", "One row per sample. raw_* columns come from the unfiltered alignment, "
                        "final_* columns from the BAM used downstream, dup_* columns from Picard MarkDuplicates. NA means the "
                        "source file does not exist (for example pre-aligned BAMs have no raw statistics)."),
    "fragment_sizes": ("Fragment sizes", "One row per sample: summary of the fragment length of properly paired reads. "
                       "pct_lt_120, pct_120_300 and pct_gt_300 are the shares of fragments below 120 bp, from 120 to 300 bp, "
                       "and above 300 bp. The size histogram is in fragment_size_distribution.tsv."),
    "library_complexity": ("Library complexity", "NRF, PBC1 and PBC2 from fragment coordinates; elc_* columns from "
                           "Picard EstimateLibraryComplexity; markdup_* columns from Picard MarkDuplicates."),
    "frip": ("Fraction of reads in peaks (FRiP)", "One row per sample and peak set. frip = reads_in_peaks / "
             "total_mapped_alignments; each sample is scored against the peak set named in peak_source."),
    "peak_stats": ("Peak statistics", "One row per peak file: number of peaks, width distribution, and the share of peaks "
                   "on numbered chromosomes."),
    "tss_enrichment": ("TSS enrichment", "enrichment = signal at the TSS (highest bin within the centre window) divided by "
                       "the mean of the outer flanks; max_over_min is the legacy score (profile maximum over minimum)."),
    "correlation_matrix": ("Correlation matrix", "Spearman coefficients behind the heatmap."),
    "replicate_correlation": ("Replicate concordance", "One row per pair of libraries of the same group; groups with one "
                              "library have status single_library and no values. pearson_all and spearman_all use all "
                              "non-blacklisted bins with reads; the *_signal columns use the top signal_frac of those bins "
                              "by mean log2 CPM. nonrep_pearson_all_median is the median pearson_all of either library "
                              "with the libraries of the same antibody in other groups (nonrep_n pairs). Descriptive "
                              "values, no thresholds."),
}

SECTIONS = [
    ("alignment", "Alignment", ["alignment_stats"], lambda s: s.startswith("alignment_")),
    ("fragments", "Fragment sizes", ["fragment_sizes"], lambda s: s.startswith("fragment_sizes_")),
    ("complexity", "Library complexity", ["library_complexity"], lambda s: s in ("library_complexity", "preseq_ccurve")),
    ("fingerprint", "Fingerprint", [], lambda s: s.startswith("fingerprint_")),
    ("correlation", "Sample correlation", ["correlation_matrix"], lambda s: s.startswith("correlation_")),
    ("repcor", "Replicate concordance", ["replicate_correlation"], lambda s: s.startswith("replicate_")),
    ("frip", "FRiP", ["frip"], lambda s: False),
    ("peaks", "Peaks", ["peak_stats"], lambda s: False),
    ("tss", "TSS enrichment", ["tss_enrichment"], lambda s: s.endswith(".tss_profile")),
    ("norm_qc", "Normalisation QC", [], lambda s: re.search(r"\.ma(\.p\d+)?$", s) is not None),
]
ORDER = {k: i for i, k in enumerate(PLOTS)}

CSS = """
:root{--bg:#ffffff;--fg:#222222;--muted:#5f6670;--rule:#dfe3e7;--accent:#86187f;--card:#f7f8fa}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15171a;--fg:#e8eaed;--muted:#a7b4c2;--rule:#333a42;--accent:#c58ac1;--card:#1e2226}}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:21px;margin:40px 0 8px;border-bottom:2px solid var(--rule);padding-bottom:4px}
h3{font-size:17px;margin:28px 0 6px}.sub{color:var(--muted);margin:0 0 16px}
nav ul{columns:2;padding-left:18px}nav a{color:var(--accent)}
.tablewrap{overflow-x:auto;margin:8px 0 16px}
table{border-collapse:collapse;font-size:13px}th,td{border:1px solid var(--rule);padding:3px 8px;text-align:right;white-space:nowrap}
th{background:var(--card)}td:first-child,th:first-child{text-align:left}
figure{margin:12px 0;padding:8px;background:#fff;border:1px solid var(--rule);border-radius:6px;overflow-x:auto}
figcaption{font-size:12px;color:#5f6670;margin-top:4px}
svg.qc-svg{max-width:100%;height:auto;display:block;margin:0 auto}
.note{background:var(--card);border-left:4px solid var(--accent);padding:8px 14px;margin:10px 0;border-radius:4px}
.note p{margin:4px 0}.note ul{margin:4px 0 4px 18px;padding:0}
.label{font-weight:600}
"""


def esc(s):
    return html.escape(str(s), quote=True)


def fmt_cell(v):
    try:
        f = float(v)
    except ValueError:
        return esc(v)
    if f == int(f) and "." not in v and "e" not in v.lower():
        return esc(v)
    return esc("%.4g" % f)


def read_tsv(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    return [l.split("\t") for l in lines]


def table_html(rows):
    if not rows:
        return "<p>The table is empty.</p>"
    out = ['<div class="tablewrap"><table><thead><tr>']
    out += ["<th>%s</th>" % esc(h) for h in rows[0]]
    out.append("</tr></thead><tbody>")
    for r in rows[1:]:
        out.append("<tr>" + "".join("<td>%s</td>" % fmt_cell(c) for c in r) + "</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


SVG_START = re.compile(r"<svg\b", re.IGNORECASE)


def namespace_svg(text, prefix):
    """Return the <svg>...</svg> element of `text` with every id renamed to '<prefix>-<id>' and every reference to
    those ids (href="#id", xlink:href, url(#id), CSS #id selectors, begin/end timing refs excluded) rewritten."""
    m = SVG_START.search(text)
    end = text.lower().rfind("</svg>")
    if not m or end < 0:
        return None
    svg = text[m.start():end + len("</svg>")]
    ids = set(re.findall(r"""\bid\s*=\s*["']([^"']+)["']""", svg))

    def new(i):
        return "%s-%s" % (prefix, i)

    svg = re.sub(r"""(\bid\s*=\s*)(["'])([^"']+)\2""", lambda mm: "%s%s%s%s" % (mm.group(1), mm.group(2), new(mm.group(3)), mm.group(2)), svg)
    svg = re.sub(r"""(href\s*=\s*)(["'])#([^"']+)\2""",
                 lambda mm: "%s%s#%s%s" % (mm.group(1), mm.group(2), new(mm.group(3)) if mm.group(3) in ids else mm.group(3), mm.group(2)), svg)
    svg = re.sub(r"""url\(\s*(["']?)#([^)"'\s]+)\1\s*\)""",
                 lambda mm: "url(%s#%s%s)" % (mm.group(1), new(mm.group(2)) if mm.group(2) in ids else mm.group(2), mm.group(1)), svg)

    def style_fix(mm):
        body = mm.group(2)
        body = re.sub(r"#([A-Za-z_][\w.\-:]*)", lambda k: "#" + (new(k.group(1)) if k.group(1) in ids else k.group(1)), body)
        return mm.group(1) + body + mm.group(3)

    svg = re.sub(r"(<style\b[^>]*>)(.*?)(</style>)", style_fix, svg, flags=re.S | re.IGNORECASE)
    # responsive sizing: keep viewBox, add class
    first_tag_end = svg.find(">")
    tag = svg[:first_tag_end]
    if "class=" in tag:
        tag = re.sub(r"""class\s*=\s*(["'])""", r'class=\1qc-svg ', tag, count=1)
    else:
        tag += ' class="qc-svg"'
    return tag + svg[first_tag_end:]


def doc_for(stem):
    if stem in PLOTS:
        return stem, PLOTS[stem]
    if stem.startswith("fingerprint_"):
        return "fingerprint", PLOTS["fingerprint"]
    if stem.endswith(".tss_profile"):
        return "tss_profile", PLOTS["tss_profile"]
    if re.search(r"\.ma(\.p\d+)?$", stem):
        return "norm_ma", PLOTS["norm_ma"]
    if stem.startswith("replicate_scatter"):
        return "replicate_correlation", PLOTS["replicate_correlation"]
    if stem.startswith("correlation_"):
        return "correlation_heatmap", PLOTS["correlation_heatmap"]
    return stem, GENERIC


def text_block(d):
    out = ['<div class="note">']
    out.append('<p><span class="label">What it measures.</span> %s</p>' % esc(d["measures"]))
    out.append('<p><span class="label">Ideal case.</span> %s</p>' % esc(d["ideal"]))
    out.append('<p><span class="label">What a deviation means.</span></p><ul>%s</ul>' % "".join("<li>%s</li>" % esc(x) for x in d["deviations"]))
    if d["edge"]:
        out.append('<p><span class="label">Edge cases.</span></p><ul>%s</ul>' % "".join("<li>%s</li>" % esc(x) for x in d["edge"]))
    out.append("</div>")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--tables-dir", default=".")
    ap.add_argument("--plots-dir", default=".")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    tables = {}
    if os.path.isdir(a.tables_dir):
        for f in sorted(os.listdir(a.tables_dir)):
            stem = f[:-4] if f.endswith(".tsv") else None
            if stem in TABLES:
                tables[stem] = os.path.join(a.tables_dir, f)
    plots = []
    if os.path.isdir(a.plots_dir):
        plots = sorted(f for f in os.listdir(a.plots_dir) if f.endswith(".svg"))
    if not tables and not plots:
        sys.stderr.write("ERROR: no QC table or SVG found in %s / %s\n" % (a.tables_dir, a.plots_dir))
        sys.exit(3)

    body, nav, counter, used = [], [], 0, set()
    for sid, title, tnames, matcher in SECTIONS:
        sec_tables = [t for t in tnames if t in tables]
        sec_plots = [p for p in plots if matcher(p[:-4]) and p not in used]
        if not sec_tables and not sec_plots:
            continue
        used.update(sec_plots)
        nav.append('<li><a href="#%s">%s</a></li>' % (sid, esc(title)))
        body.append('<section id="%s"><h2>%s</h2>' % (sid, esc(title)))
        for t in sec_tables:
            name, blurb = TABLES[t]
            body.append("<h3>%s</h3><p>%s</p>%s" % (esc(name), esc(blurb), table_html(read_tsv(tables[t]))))
        # plots grouped by text kind, one text block per kind, kinds in registry order
        kinds, docs = {}, {}
        for p in sec_plots:
            key, d = doc_for(p[:-4])
            kinds.setdefault(key, []).append(p)
            docs[key] = d
        for key, files in sorted(kinds.items(), key=lambda kv: (ORDER.get(kv[0], 999), kv[0])):
            d = docs[key]
            body.append("<h3>%s</h3>%s" % (esc(d["title"]), text_block(d)))
            for p in files:
                with open(os.path.join(a.plots_dir, p), encoding="utf-8", errors="replace") as fh:
                    raw = fh.read()
                counter += 1
                svg = namespace_svg(raw, "p%02d" % counter)
                if svg is None:
                    body.append("<p>%s could not be embedded (no &lt;svg&gt; element).</p>" % esc(p))
                    continue
                body.append("<figure>%s<figcaption>%s</figcaption></figure>" % (svg, esc(p)))
        body.append("</section>")

    leftovers = [p for p in plots if p not in used]
    if leftovers:
        nav.append('<li><a href="#other">Other plots</a></li>')
        body.append('<section id="other"><h2>Other plots</h2>')
        for p in leftovers:
            with open(os.path.join(a.plots_dir, p), encoding="utf-8", errors="replace") as fh:
                raw = fh.read()
            counter += 1
            svg = namespace_svg(raw, "p%02d" % counter)
            if svg:
                key, d = doc_for(p[:-4])
                body.append("<h3>%s</h3>%s<figure>%s<figcaption>%s</figcaption></figure>" % (esc(p), text_block(d), svg, esc(p)))
        body.append("</section>")

    doc = ("<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
           "<title>%s QC report</title><style>%s</style></head><body><main>"
           "<h1>%s QC report</h1><p class=\"sub\">Alignment, fragment, complexity, enrichment and peak quality control.</p>"
           "<nav><ul>%s</ul></nav>%s</main></body></html>\n") % (esc(a.project), CSS, esc(a.project), "".join(nav), "".join(body))
    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, "%s_qc_report.html" % a.project)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print("wrote %s (%d tables, %d plots)" % (out, len(tables), len(plots)))


if __name__ == "__main__":
    main()
