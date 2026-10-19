#!/usr/bin/env python3
"""
discover_anchor.py
=====================
Finds the constant vector sequence between a sample's barcode and its guide
(sgRNA) — the "anchor" guide_anchor_trim.sh then searches for, per read, to
position each read's guide exactly regardless of any variable-length
spacing before or after the anchor. Requires no prior knowledge of the
anchor's sequence OR length: it's reconstructed directly from the reads.

Method: greedy consensus extension. Demultiplexed reads start right after
the barcode (true for any plasmid design — that's what demultiplexing does,
not something specific to one library), so position 0 is where a constant
region should begin if this amplicon has one. For each length k = 4, 5, 6,
..., find the most common exact k-mer at position [0, k) across a sample of
reads, and track what fraction of reads carry it. A true constant base is
retained in ~95-99%+ of reads that had the previous k-1 bases right (lost
only to ordinary sequencing error); once k walks past the constant region
into the guide (itself not perfectly random — e.g. many sgRNA designs are
enriched for a 5' G for U6 promoter efficiency), the retention ratio
between consecutive k collapses sharply. That collapse point is the
anchor's end.

Library-aware refinement (--library): consensus extension alone can't tell
a constant vector base from a guide base that happens to be the same in
(nearly) every guide — e.g. a library whose 20-mers all start with the
U6 "G" (G + 19 nt): every read then carries that G right after the vector,
so the consensus anchor swallows the guide's first base, and every trimmed
read is shifted by one (nothing maps). So, given this library's guides, the
consensus anchor is trimmed back 0, 1, 2, ... --max-refine-trim bases from
its 3' end, and the candidate kept is the one followed by an exact library
guide in the most sampled reads (ties -> the longer candidate). If no
candidate is followed by any library guide (e.g. this sample holds none of
this library's reads), the consensus anchor is kept unchanged and that is
recorded in the refinement table.

USAGE:
  python3 discover_anchor.py --fastq=<sample.fastq.gz> --project-name=<name> \
      --library-name=<name> [--library=<lib_mageck_library.txt>] \
      [--n-reads=200000] [--min-length=8] [--max-length=40] \
      [--min-frequency=30] [--ratio-cutoff=0.70] [--max-refine-trim=8] \
      [--outdir=<dir>]

  --library   MAGeCK-format library (tab-separated, header, sequence in
              column 2 — build_mageck_library.py's output). Optional;
              without it no refinement is done and the refinement table
              records that.

OUTPUTS (<prefix> = <project-name>_<library-name>):
  <prefix>_anchor.txt               the discovered anchor sequence, one line
  <prefix>_anchor_discovery.csv     k, top k-mer, count, % of reads, ratio vs k-1
                                    — the full discovery trace, for QC/troubleshooting
  <prefix>_anchor_refinement.csv    one row per candidate tried: bases trimmed,
                                    candidate, reads starting with it, reads with an
                                    exact library guide right after it, chosen (TRUE/FALSE)

Exits 3 (a data problem; the pipeline doesn't retry it) if no anchor at
least --min-length long can be discovered with at least --min-frequency%
of reads agreeing on it — e.g. the sample truly has no constant region
here, or reads are lower quality than assumed.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import Counter


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(3)


def open_fastq(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def sample_reads(path: str, n_reads: int) -> list[str]:
    reads = []
    with open_fastq(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            if line_no % 4 != 2:
                continue
            reads.append(line.strip())
            if len(reads) >= n_reads:
                break
    return reads


def discover(reads: list[str], min_length: int, max_length: int, ratio_cutoff: float):
    """Returns (anchor, anchor_freq, trace_rows). trace_rows is a list of
    dicts for the QC CSV; anchor is None if no candidate reached
    min_length. Stops (and falls back to the last accepted candidate)
    whenever the retention ratio collapses, max_length is reached, or no
    read is long enough to extend further — the three ways this loop can
    end all need the same "use the last good candidate" handling."""
    trace = []
    prev_freq, prev_seq = None, None
    n = len(reads)
    for k in range(4, max_length + 1):
        counts = Counter(r[:k] for r in reads if len(r) >= k)
        if not counts:
            break
        top_seq, top_count = counts.most_common(1)[0]
        freq = top_count / n * 100 if n else 0.0
        ratio = (freq / prev_freq) if prev_freq else 1.0
        trace.append({"k": k, "top_kmer": top_seq, "count": top_count,
                       "pct_of_reads": round(freq, 3), "ratio_vs_prev": round(ratio, 4)})
        if prev_freq is not None and ratio < ratio_cutoff:
            break
        prev_freq, prev_seq = freq, top_seq

    if prev_seq and len(prev_seq) >= min_length:
        return prev_seq, prev_freq, trace
    return None, None, trace


def read_library(path: str) -> set[str]:
    """Guide sequences from a MAGeCK-format library (header line; sequence = column 2)."""
    guides = set()
    with open(path, encoding="utf-8-sig") as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) >= 2 and parts[1].strip():
                guides.add(parts[1].strip().upper())
    return guides


def refine(reads: list[str], anchor: str, guides: set[str], max_trim: int, floor: int):
    """Score anchor[:len-t] for t = 0..max_trim by how many reads start with it AND
    carry an exact library guide immediately after it. Returns (best_anchor, rows)."""
    lengths = sorted({len(g) for g in guides}, reverse=True)
    rows = []
    for t in range(0, max_trim + 1):
        cand = anchor[: len(anchor) - t]
        if len(cand) < floor:
            break
        n = len(cand)
        with_prefix = hits = 0
        for r in reads:
            if not r.startswith(cand):
                continue
            with_prefix += 1
            if any(r[n:n + L] in guides for L in lengths if len(r) >= n + L):
                hits += 1
        rows.append({"bases_trimmed": t, "candidate": cand, "length": n,
                     "reads_starting_with_candidate": with_prefix,
                     "reads_with_library_guide_after": hits,
                     "pct_of_sampled_reads": round(hits / len(reads) * 100, 3) if reads else 0.0,
                     "chosen": False})
    best = max(rows, key=lambda x: (x["reads_with_library_guide_after"], x["length"])) if rows else None
    if best is None or best["reads_with_library_guide_after"] == 0:
        for x in rows:
            x["chosen"] = x["bases_trimmed"] == 0
        return anchor, rows, False
    best["chosen"] = True
    return best["candidate"], rows, True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--project-name", required=True)
    ap.add_argument("--library-name", required=True)
    ap.add_argument("--n-reads", type=int, default=200_000)
    ap.add_argument("--min-length", type=int, default=8,
                     help="shorter than this isn't a specific enough anchor to trust (default 8)")
    ap.add_argument("--max-length", type=int, default=40,
                     help="how far into the read to search (default 40)")
    ap.add_argument("--min-frequency", type=float, default=30.0,
                     help="minimum %% of reads that must agree on the final anchor (default 30)")
    ap.add_argument("--ratio-cutoff", type=float, default=0.70,
                     help="stop extending when retention drops below this fraction of the previous length's rate (default 0.70)")
    ap.add_argument("--library", default=None,
                     help="MAGeCK-format library of this library's guides; enables library-aware refinement")
    ap.add_argument("--max-refine-trim", type=int, default=8,
                     help="how many 3' bases refinement may trim off the consensus anchor (default 8)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    prefix = f"{args.project_name}_{args.library_name}"

    reads = sample_reads(args.fastq, args.n_reads)
    if not reads:
        fail(f"{args.fastq} has no reads to sample.")
    print(f"Sampled {len(reads):,} reads from {args.fastq}", file=sys.stderr)

    anchor, anchor_freq, trace = discover(reads, args.min_length, args.max_length, args.ratio_cutoff)

    os.makedirs(args.outdir, exist_ok=True)
    trace_path = os.path.join(args.outdir, f"{prefix}_anchor_discovery.csv")
    with open(trace_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["k", "top_kmer", "count", "pct_of_reads", "ratio_vs_prev"])
        writer.writeheader()
        writer.writerows(trace)
    print(f"  Saved: {trace_path}", file=sys.stderr)

    if anchor is None:
        fail(
            f"could not discover a constant anchor at least {args.min_length}bp long in "
            f"{args.fastq} (sampled {len(reads):,} reads). This usually means either this "
            f"library's amplicon has no constant region right after the barcode, or "
            f"{args.fastq} isn't actually demultiplexed/positioned as expected. See "
            f"{trace_path} for the full retention trace. If this library genuinely has no "
            f"anchor, set trim5 in config.yaml to bypass discovery and use a fixed offset instead."
        )
    if anchor_freq < args.min_frequency:
        fail(
            f"the best candidate anchor ({anchor}, length {len(anchor)}) was only supported "
            f"by {anchor_freq:.1f}% of sampled reads, below the {args.min_frequency}% "
            f"confidence floor. See {trace_path} for the full retention trace — this usually "
            f"means the reads are noisier than expected, or {args.fastq} isn't what this "
            f"library's samples should look like."
        )

    refine_path = os.path.join(args.outdir, f"{prefix}_anchor_refinement.csv")
    fields = ["bases_trimmed", "candidate", "length", "reads_starting_with_candidate",
              "reads_with_library_guide_after", "pct_of_sampled_reads", "chosen"]
    if args.library:
        guides = read_library(args.library)
        if not guides:
            fail(f"--library {args.library} holds no guide sequences.")
        consensus = anchor
        anchor, refine_rows, used = refine(reads, consensus, guides, args.max_refine_trim,
                                           min(args.min_length, len(consensus)))
        if not used:
            print(f"  [refinement] no sampled read has an exact {args.library_name} guide right after "
                  f"any candidate — keeping the consensus anchor {consensus}", file=sys.stderr)
        elif anchor != consensus:
            print(f"  [refinement] consensus anchor {consensus} (length {len(consensus)}) ends inside the "
                  f"guide: trimmed to {anchor} (length {len(anchor)}), which is followed by a library "
                  f"guide in the most reads", file=sys.stderr)
    else:
        refine_rows = [{"bases_trimmed": 0, "candidate": anchor, "length": len(anchor),
                        "reads_starting_with_candidate": "", "reads_with_library_guide_after": "",
                        "pct_of_sampled_reads": "", "chosen": "TRUE (no --library: refinement not run)"}]
    with open(refine_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(refine_rows)
    print(f"  Saved: {refine_path}", file=sys.stderr)

    anchor_path = os.path.join(args.outdir, f"{prefix}_anchor.txt")
    with open(anchor_path, "w", encoding="utf-8") as fh:
        fh.write(anchor + "\n")
    print(f"  Saved: {anchor_path}", file=sys.stderr)
    print(f"\nDiscovered anchor for {args.library_name}: {anchor} "
          f"(length {len(anchor)}, {anchor_freq:.1f}% of sampled reads)", file=sys.stderr)


if __name__ == "__main__":
    main()
