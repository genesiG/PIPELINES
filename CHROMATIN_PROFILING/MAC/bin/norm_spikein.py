#!/usr/bin/env python3
"""norm_spikein.py -- physical spike-in (E. coli) scale factors for ONE normalisation run.

CLI contract (called by the NORM_SPIKEIN module; runs in cutandrun-py, which has no samtools/pysam, so the read counts
come from `samtools idxstats` text files produced upstream):

  norm_spikein.py --target-idxstats T1.idxstats T2.idxstats ...   idxstats of the FINAL target (human) BAMs
                  --spike-idxstats  S1.idxstats S2.idxstats ...   idxstats of the SPIKE-IN BAMs (deduplicated when
                                                                  spikein_dedup is true; that is decided upstream)
                  --sample-ids      id1 id2 ...                   same order and length as the two idxstats lists
                  --norm-run        <name>                        run label (goes into the note column)
                  --out             <name>.scale_factors.tsv
                  [--spike-prefix   ecoli_chr]                    optional: count only spike contigs in the spike idxstats
                                                                  and exclude them from the target idxstats (default:
                                                                  "" = sum every contig of each file)
                  [--min-spike-reads 1000]                        optional: qc_pass = FALSE below this many spike reads
                  [--fail-soft]                                   optional flag (diagnostic use, spike-in is not the run's
                                                                  primary method): a library with 0 spike-in reads does not
                                                                  stop the script; every row is written with NA factors,
                                                                  qc_pass FALSE and note "UNAVAILABLE: <id> has 0 spike-in
                                                                  reads" (ids joined by "; "), exit 0. Without the flag
                                                                  (primary use) S = 0 exits 3 as before. A target library with
                                                                  0 reads is still a data problem (exit 3).
Every list option accepts space-separated values and/or comma-separated values.

Output <name>.scale_factors.tsv: sample_id, method, raw_factor, bigwig_multiplier, edger_norm_factor, qc_pass, note
Exit codes: 0 ok (also the UNAVAILABLE table of --fail-soft); 2 usage error; 3 data problem (unreadable idxstats, zero target
or zero spike reads without --fail-soft, length mismatch).

Statistic and assumptions
  Question: how much more/less antibody-targeted material does each library hold per unit of spike-in chromatin?
  Assumption: the same amount of spike-in chromatin was added per cell/sample, the spike-in is captured with the same
  efficiency in every library, and the spike-in reads in the (deduplicated) spike BAM are not saturated. It does NOT require
  most genomic regions to be unchanged, so it is the right tool when a global change is expected (e.g. EZH2/EED inhibitor).
  Low spike-in counts make the factor noisy (relative SE ~ 1/sqrt(S)), hence the --min-spike-reads flag.

Algebra (unchanged from legacy step_4c, EpiCypher CUTANA protocol)
  S = spike-in reads (sum of column 3 of the spike idxstats), T = target reads (sum of column 3 of the target idxstats)
  raw_factor         = (S + T) / S = 1 / spike_fraction
  bigwig_multiplier  = raw_factor / min(raw_factor over the run)        (legacy bamCov_SF; >= 1, anchored at the library with
                       the highest spike fraction). bamCoverage multiplies --scaleFactor into the CPM/RPKM value, and RPM
                       already divides by T, so RPM x raw ~ (T/S) x counts/T = counts/S: the signal is expressed per
                       spike-in read. Use with --normalizeUsing CPM or RPKM.
  edger_norm_factor  = (min(raw)/raw) / geometric_mean(min(raw)/raw)     (edgeR convention: norm.factors have geometric
                       mean 1; effective library size = lib.size x norm.factors, so CPM = count/(lib.size x nf)).
                       min(raw)/raw is the legacy "SF" column; it is proportional to S/(S+T), i.e. to S/T up to a factor
                       (1 + S/T) that is below 0.3 percent for the NUC libraries (raw factors 390-1160, S/T = 1/(raw-1)).
  Invariant: bigwig_multiplier x edger_norm_factor = 1 / geometric_mean(min(raw)/raw), the same for every library of the run.
  Library sizes: edgeR's lib.size must be the target-read total in the same BAM universe (final target BAM) for nf to mean
  what is stated; csaw `totals` of the same BAMs are within a small fraction of the idxstats mapped-read sum.
"""

import argparse
import math
import os
import sys

COLUMNS = ["sample_id", "method", "raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "note"]


def die(msg, status=3):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(status)


def expand(values):
    out = []
    for v in values:
        out.extend(x.strip() for x in v.split(",") if x.strip())
    return out


def read_idxstats(path, spike_prefix, want_spike):
    """Sum column 3 (mapped reads). With a prefix: spike file -> contigs starting with it; target file -> all others."""
    if not os.path.isfile(path):
        die("idxstats file not found: %s" % path)
    total = 0
    n_lines = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            n_lines += 1
            if f[0] == "*":
                continue
            try:
                mapped = int(f[2])
            except ValueError:
                die("%s: non-integer mapped count in line '%s'" % (path, line.strip()))
            if spike_prefix:
                is_spike = f[0].startswith(spike_prefix)
                if is_spike != want_spike:
                    continue
            total += mapped
    if n_lines == 0:
        die("%s is empty or not a samtools idxstats file" % path)
    return total


def main():
    ap = argparse.ArgumentParser(description="Spike-in scale factors (raw=(S+T)/S, min-anchored) from idxstats files.")
    ap.add_argument("--target-idxstats", nargs="+", required=True)
    ap.add_argument("--spike-idxstats", nargs="+", required=True)
    ap.add_argument("--sample-ids", nargs="+", required=True)
    ap.add_argument("--norm-run", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--spike-prefix", default="")
    ap.add_argument("--min-spike-reads", type=int, default=1000)
    ap.add_argument("--fail-soft", action="store_true")
    a = ap.parse_args()

    target_files = expand(a.target_idxstats)
    spike_files = expand(a.spike_idxstats)
    ids = expand(a.sample_ids)
    if not (len(target_files) == len(spike_files) == len(ids)):
        die("list lengths differ: %d target idxstats, %d spike idxstats, %d sample ids"
            % (len(target_files), len(spike_files), len(ids)), 2)
    if len(ids) < 2:
        die("a spike-in normalisation run needs at least 2 libraries (got %d)" % len(ids), 2)
    if len(set(ids)) != len(ids):
        die("duplicated sample ids: %s" % ", ".join(sorted({i for i in ids if ids.count(i) > 1})), 2)

    S, T, problems, zero_spike = [], [], [], []
    for sid, tf, sf in zip(ids, target_files, spike_files):
        for path in (tf, sf):  # soft check that the file belongs to the sample (names are staged as <sample>.*)
            if sid not in os.path.basename(path):
                sys.stderr.write("WARNING: sample id '%s' not found in file name '%s'; check the pairing order\n"
                                 % (sid, os.path.basename(path)))
        t = read_idxstats(tf, a.spike_prefix, want_spike=False)
        s = read_idxstats(sf, a.spike_prefix, want_spike=True)
        if t == 0:
            problems.append("%s: 0 mapped reads in the target idxstats" % sid)
        if s == 0:
            zero_spike.append(sid)
            problems.append("%s: 0 mapped reads in the spike-in idxstats (spike fraction 0 -> infinite factor)" % sid)
        S.append(s)
        T.append(t)
    if zero_spike and a.fail_soft and not [p for p in problems if "0 mapped reads in the target" in p]:
        reason = "; ".join("%s has 0 spike-in reads" % sid for sid in zero_spike)
        sys.stderr.write("WARNING: spike-in unavailable (--fail-soft): %s\n" % reason)
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write("\t".join(COLUMNS) + "\n")
            for sid in ids:
                fh.write("\t".join([sid, "spikein", "NA", "NA", "NA", "FALSE",
                                    "UNAVAILABLE: %s; run=%s" % (reason, a.norm_run)]) + "\n")
        print("Wrote %s (UNAVAILABLE: %s)" % (a.out, reason))
        return
    if problems:
        die("cannot compute spike-in factors:\n  " + "\n  ".join(problems))

    raw = [(s + t) / s for s, t in zip(S, T)]
    raw_min = min(raw)
    anchor = ids[raw.index(raw_min)]
    bigwig = [r / raw_min for r in raw]
    nf0 = [raw_min / r for r in raw]
    gm = math.exp(sum(math.log(x) for x in nf0) / len(nf0))
    nf = [x / gm for x in nf0]

    low = []
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sid, s, t, r, b, n in zip(ids, S, T, raw, bigwig, nf):
            ok = s >= a.min_spike_reads
            if not ok:
                low.append(sid)
            note = ("run=%s; spike_reads=%d; target_reads=%d; spike_fraction=%.6g; anchor=%s%s"
                    % (a.norm_run, s, t, s / (s + t), anchor,
                       "" if ok else "; LOW SPIKE-IN COUNT (< %d): factor noisy" % a.min_spike_reads))
            fh.write("\t".join([sid, "spikein", repr(r), repr(b), repr(n), "TRUE" if ok else "FALSE", note]) + "\n")
    if low:
        sys.stderr.write("WARNING: spike-in reads below --min-spike-reads (qc_pass=FALSE): %s\n" % ", ".join(low))
    print("Wrote %s (%d libraries; anchor = %s; min raw factor %.6f)" % (a.out, len(ids), anchor, raw_min))


if __name__ == "__main__":
    main()
