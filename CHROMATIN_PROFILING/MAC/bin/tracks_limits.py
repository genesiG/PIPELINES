#!/usr/bin/env python3
r"""tracks_limits.py -- shared y maxima of the grouped bigWig tracks of one locus, written into a resolved ini.

Called by process TRACKS_LIMITS (one call per track set and locus):
    tracks_limits.py --ini NAME.tracks.ini --region chr1:1000-2000 --out-prefix SLUG [--limit-stat mean_bins]

Arguments:
  --ini         tracks.ini from TRACKS_INI (relative bigWig file names, resolved from the working directory)
  --region      chr:start-end of the plotted locus; commas inside the numbers are allowed
  --out-prefix  SLUG (the track set and locus slug); outputs SLUG.resolved.ini and SLUG_limits.tsv
  --limit-stat  mean_bins (default) | max_bin: the statistic of each track that the group limit is based on
                (params.tracks_limit_stat)

Algorithm: every ini section whose max_value is 'group:NAME' gives its file= and number_of_bins= (default 700).
Per track, the track peak depends on --limit-stat:
  mean_bins  bw.stats(chrom, start, end, nBins=bins, type="mean") over the region gives the per-bin means
             (the end is clipped to the chromosome length; None bins are ignored). Track peak = the largest bin mean
             (0 if none). EED rule.
  max_bin    v = bw.stats(chrom, start, end, type="max")[0], the maximum over the whole locus window (the end clipped
             as above; the bigWig's own bins). Track peak = v (0 if v is None or NaN). CBX2 rule.
Per group: peak = the largest track peak of its tracks (0 if none). Limit = round_up_2sf(1.05 * peak) in both modes;
round_up_2sf rounds up to 2 significant figures (a non-positive input gives 1).
pyBigWig is imported only when the ini has at least one grouped track.

Outputs (working directory):
  SLUG.resolved.ini  the ini with every 'max_value = group:NAME' line replaced by 'max_value = <limit>'; all other
                     lines are unchanged. With no groups it is an exact copy of the ini.
  SLUG_limits.tsv    header 'group<TAB>track<TAB>file<TAB>stat<TAB>track_peak<TAB>limit', one row per grouped track
                     (track = the ini section name; stat = the --limit-stat mode; track_peak = the track's value, %.6g).
                     Header only when there are no groups.
Stdout: none. Stderr: one line per group with its shared maximum and the statistic.
Exit 3: the ini or a bigWig is missing or not a bigWig; a chromosome of the region is not in a grouped bigWig;
        the region is malformed or lies outside its chromosome; a 'group:' max_value remains in the resolved ini;
        a grouped section has no file=.
Exit 2: a required argument is missing.
"""
import argparse
import math
import os
import re
import sys

DEFAULT_BINS = 700
REGION_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")
KEY_RE = re.compile(r"^\s*([A-Za-z_]+)\s*=\s*(\S*)")
GROUP_LINE_RE = re.compile(r"^\s*max_value\s*=\s*group:")


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def round_up_2sf(x):
    """Round x up to 2 significant figures (1 if x is not positive)."""
    if not x > 0:
        return 1
    e = math.floor(math.log10(x)) - 1
    return round(math.ceil(round(x / 10**e, 9)) * 10**e, 10)


def parse_sections(lines):
    """Return the list of ini sections as dicts: name, file, bins, group, line (index of the max_value line)."""
    sections, cur = [], None
    for i, line in enumerate(lines):
        if line.startswith("["):
            cur = {"name": line.strip()[1:-1], "file": None, "bins": DEFAULT_BINS, "group": None, "line": None}
            sections.append(cur)
            continue
        if cur is None:
            continue
        m = KEY_RE.match(line)
        if not m:
            continue
        key, val = m.group(1), m.group(2)
        if key == "file":
            cur["file"] = val
        elif key == "number_of_bins":
            try:
                cur["bins"] = int(val)
            except ValueError:
                die("section [%s]: number_of_bins is not an integer ('%s')" % (cur["name"], val))
        elif key == "max_value" and val.startswith("group:"):
            cur["group"] = val[len("group:"):]
            cur["line"] = i
    return sections


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ini", required=True)
    ap.add_argument("--region", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--limit-stat", choices=["mean_bins", "max_bin"], default="mean_bins")
    a = ap.parse_args()

    try:
        with open(a.ini, encoding="utf-8") as fh:
            text = fh.read()
    except OSError as e:
        die("cannot read ini %s: %s" % (a.ini, e))
    lines = text.split("\n")
    out_ini = a.out_prefix + ".resolved.ini"
    out_tsv = a.out_prefix + "_limits.tsv"

    sections = parse_sections(lines)
    grouped = [sec for sec in sections if sec["group"] is not None]
    rows = []
    limit = {}

    if grouped:
        for sec in grouped:
            if not sec["file"]:
                die("section [%s] has max_value = group:%s but no file" % (sec["name"], sec["group"]))
            if not os.path.isfile(sec["file"]):
                die("bigWig not found: %s (section [%s])" % (sec["file"], sec["name"]))

        region = a.region.replace(",", "")
        m = REGION_RE.match(region)
        if not m:
            die("region is not chr:start-end (got '%s')" % a.region)
        chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

        import pyBigWig

        peak = {}
        for sec in grouped:
            try:
                bw = pyBigWig.open(sec["file"])
            except RuntimeError:
                bw = None
            if bw is None:
                die("cannot open bigWig %s (section [%s])" % (sec["file"], sec["name"]))
            try:
                chroms = bw.chroms()
                if chrom not in chroms:
                    die("chromosome %s is not in bigWig %s (section [%s])" % (chrom, sec["file"], sec["name"]))
                clip_end = min(end, int(chroms[chrom]))
                if start >= clip_end:
                    die("region %s lies outside chromosome %s (length %d)" % (a.region, chrom, chroms[chrom]))
                if a.limit_stat == "max_bin":
                    v = bw.stats(chrom, start, clip_end, type="max")[0]
                    track_peak = 0.0 if v is None or math.isnan(v) else v
                else:
                    vals = bw.stats(chrom, start, clip_end, nBins=sec["bins"], type="mean")
                    finite = [v for v in vals if v is not None and not math.isnan(v)]
                    track_peak = max(finite) if finite else 0.0
            finally:
                bw.close()
            sec["track_peak"] = track_peak
            peak[sec["group"]] = max(peak.get(sec["group"], 0), track_peak)

        limit = {g: round_up_2sf(1.05 * v) for g, v in peak.items()}
        for sec in grouped:
            lines[sec["line"]] = "max_value = %g" % limit[sec["group"]]
            rows.append((sec["group"], sec["name"], sec["file"], a.limit_stat, sec["track_peak"], limit[sec["group"]]))

    resolved = "\n".join(lines)
    for line in resolved.split("\n"):
        if GROUP_LINE_RE.match(line):
            die("a 'group:' max_value remains in the resolved ini: %s" % line.strip())
    with open(out_ini, "w", encoding="utf-8") as fh:
        fh.write(resolved)
    with open(out_tsv, "w", encoding="utf-8") as fh:
        fh.write("group\ttrack\tfile\tstat\ttrack_peak\tlimit\n")
        for g, name, f, stat, mx, lim in rows:
            fh.write("%s\t%s\t%s\t%s\t%.6g\t%g\n" % (g, name, f, stat, mx, lim))

    for g in sorted(limit):
        sys.stderr.write("Shared y max %s in %s (%s): %g\n" % (g, a.region, a.limit_stat, limit[g]))


if __name__ == "__main__":
    main()
