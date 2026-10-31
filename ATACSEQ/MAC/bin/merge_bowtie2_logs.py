#!/usr/bin/env python3
"""merge_bowtie2_logs.py -- sum the paired-end bowtie2 summaries of the chunk logs of one sample.

Called by ALIGN_MERGE through bin/align_merge.sh:
    merge_bowtie2_logs.py --out S.bowtie2.log S.chunk000.bowtie2.log S.chunk001.bowtie2.log

Arguments:
  --out   output log, written in bowtie2's own summary format (only the summary block)
  LOG     one or more bowtie2 logs (bowtie2 stderr), in chunk order
Only the summary block of each log is read: from the line 'N reads; of these:' to the line 'X% overall alignment
  rate'. Every other line ([WARNING] lines, the -t timing lines) is ignored.
Counts: the counts of all logs are summed line by line. Percentages are recomputed from the summed counts as
  100 * count / parent count, 2 decimals (0.00 when the parent count is 0); the parent of a line is the nearest
  earlier line with a smaller indent that has a count. The overall rate is recomputed too: total mates = 2 * the
  first count; unaligned mates = the 'aligned 0 times' count under 'mates make up the pairs' when that section
  exists, else 2 * the 'aligned concordantly 0 times' count; rate = 100 * (total - unaligned) / total.
  A single log gives back its own summary block unchanged.
Output: the block as bowtie2 prints it ('COUNT (PCT%) text', same indentation, separators as in the input).
Exit 2: usage error (argparse). Exit 3: a log cannot be read, has no summary block, has no overall rate line,
  has unpaired reads ('were unpaired'), has a line in the block that is neither a count line nor a separator,
  or whose block differs in structure from the first log.
"""
import argparse
import re
import sys

TOP = re.compile(r"^[0-9]+ reads; of these:$")
END = re.compile(r"^[0-9.]+% overall alignment rate$")
COUNT = re.compile(r"^([ \t]*)([0-9]+)(?: \(([0-9]+\.[0-9]+)%\))? (.+)$")
SEP = re.compile(r"^([ \t]*)-+$")
MATES = "mates make up the pairs; of these:"
ZERO_MATES = "aligned 0 times"
ZERO_CONC = "aligned concordantly 0 times"


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def read_block(path):
    """Summary block of one log as a list of records: sep, indent, count, pct (bool), text, raw."""
    try:
        with open(path, encoding="utf-8") as fh:
            lines = [ln.rstrip("\r\n") for ln in fh]
    except (OSError, UnicodeDecodeError) as e:
        die("cannot read %s: %s" % (path, e))
    if any("were unpaired" in ln for ln in lines):
        die("%s has unpaired reads ('were unpaired'); only paired-end bowtie2 logs can be merged" % path)
    start = next((i for i, ln in enumerate(lines) if TOP.match(ln)), None)
    if start is None:
        die("%s: no bowtie2 summary block (no line 'N reads; of these:')" % path)
    end = next((i for i in range(start + 1, len(lines)) if END.match(lines[i])), None)
    if end is None:
        die("%s: no 'X%% overall alignment rate' line after the summary block" % path)

    recs = []
    for n, ln in enumerate(lines[start:end], start + 1):
        m = SEP.match(ln)
        if m:
            recs.append({"sep": True, "indent": m.group(1), "count": 0, "pct": False,
                         "text": ln.strip(), "raw": ln})
            continue
        m = COUNT.match(ln)
        if not m:
            die("%s line %d is neither a count line nor a separator: %s" % (path, n, ln))
        recs.append({"sep": False, "indent": m.group(1), "count": int(m.group(2)), "pct": m.group(3) is not None,
                     "text": m.group(4), "raw": ln})
    return recs


def shape(recs):
    return [(r["sep"], r["indent"], r["text"], r["pct"]) for r in recs]


def parent_of(recs, i):
    """Nearest earlier count line with a smaller indent; None when there is none."""
    for j in range(i - 1, -1, -1):
        if not recs[j]["sep"] and len(recs[j]["indent"]) < len(recs[i]["indent"]):
            return j
    return None


def find_text(recs, text, after=-1):
    return next((j for j in range(after + 1, len(recs)) if not recs[j]["sep"] and recs[j]["text"] == text), None)


def render(recs, counts, path):
    out = []
    for i, r in enumerate(recs):
        if r["sep"]:
            out.append(r["raw"])
        elif r["pct"]:
            p = parent_of(recs, i)
            if p is None:
                die("%s: line '%s' has a percentage but no parent line" % (path, r["raw"]))
            pc = counts[p]
            val = "0.00" if pc == 0 else "%.2f" % (100.0 * counts[i] / pc)
            out.append("%s%d (%s%%) %s" % (r["indent"], counts[i], val, r["text"]))
        else:
            out.append("%s%d %s" % (r["indent"], counts[i], r["text"]))
    return out


def overall_line(recs, counts):
    total = 2 * counts[0]
    m = find_text(recs, MATES)
    if m is not None:
        z = find_text(recs, ZERO_MATES, after=m)
        if z is None:
            die("the mates section has no 'aligned 0 times' line")
        unaligned = counts[z]
    else:
        z = find_text(recs, ZERO_CONC)
        if z is None:
            die("no 'aligned concordantly 0 times' line in the summary")
        unaligned = 2 * counts[z]
    rate = 100.0 * (total - unaligned) / total if total else 0.0
    return "%.2f%% overall alignment rate" % rate


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("logs", nargs="+", metavar="LOG")
    a = ap.parse_args()

    first = read_block(a.logs[0])
    counts = [0] * len(first)
    for path in a.logs:
        recs = read_block(path)
        if shape(recs) != shape(first):
            n = next((k for k, (x, y) in enumerate(zip(shape(recs), shape(first))) if x != y), min(len(recs), len(first)))
            die("%s differs in structure from %s at block line %d (%d vs %d lines in the block)"
                % (path, a.logs[0], n + 1, len(recs), len(first)))
        for i, r in enumerate(recs):
            if not r["sep"]:
                counts[i] += r["count"]

    lines = render(first, counts, a.logs[0]) + [overall_line(first, counts)]
    try:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
    except OSError as e:
        die("cannot write %s: %s" % (a.out, e))
    print("wrote %s (%d logs)" % (a.out, len(a.logs)))


if __name__ == "__main__":
    main()
