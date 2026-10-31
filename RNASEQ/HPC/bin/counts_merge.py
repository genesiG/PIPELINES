#!/usr/bin/env python3
"""counts_merge.py -- merge the per-sample featureCounts gene counts into one matrix, plus gene lengths and summaries.

Called by COUNTS_MERGE (modules/local/counts_merge.nf):
    counts_merge.py --dir . --samples a,b,c --project-name P --out-counts P_counts.txt

Arguments:
  --dir           directory with the per-sample files (all staged flat by the module):
                  <s>_counts.txt (FEATURECOUNTS output: Geneid, count; no header),
                  <s>_featurecounts.txt (featureCounts table, 2 comment/header lines first),
                  <s>_featurecounts.txt.summary (featureCounts summary)
  --samples       comma-separated sample ids in the order of the samplesheet (column order of the matrix)
  --project-name  project name (used in messages)
  --out-counts    the count matrix, e.g. P_counts.txt
Outputs (gene_length.txt and featurecounts_summary.tsv are written to the directory of --out-counts):
  P_counts.txt             header "Gene" + sample ids in --samples order; rows in featureCounts order;
                           integer counts; tab separated (the legacy step_5 layout; the legacy column
                           names carried the read suffix, here they are the sample ids)
  gene_length.txt          columns 1, 5, 6 (Geneid, Strand, Length) of the FIRST sample's featureCounts table,
                           without header (legacy: tail -n+3 | cut -f1,5,6 in every job; one copy is written here)
  featurecounts_summary.tsv  header "Status" + sample ids; one row per featureCounts summary status (Assigned,
                           Unassigned_*), the values of each sample; the summary header cell with the BAM path is
                           not used: samples are matched by file name
Checks: every sample's Geneid column must equal the first sample's (same ids, same order); the first differing
  row is named. Exit 3 for missing or empty files, a line with the wrong number of columns, a non-integer count,
  duplicated sample ids, a summary whose status rows differ from the first sample's, or a gene mismatch.
The legacy step_5 sorted the sample columns; here they follow --samples (the samplesheet order).
The merge is a plain join of identical gene columns, so the legacy inner merge is equivalent here.
"""
import argparse
import os
import sys


def die(msg):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(3)


def read_lines(path):
    if not os.path.isfile(path):
        die("missing file %s" % path)
    if os.path.getsize(path) == 0:
        die("empty file %s" % path)
    with open(path, encoding="utf-8") as fh:
        return [line.rstrip("\r\n") for line in fh]


def read_counts(path):
    genes, counts = [], []
    for n, line in enumerate(read_lines(path), 1):
        fields = line.split("\t")
        if len(fields) != 2:
            die("%s line %d: expected 2 tab-separated columns, found %d" % (path, n, len(fields)))
        try:
            counts.append(int(fields[1]))
        except ValueError:
            die("%s line %d: count '%s' is not an integer" % (path, n, fields[1]))
        genes.append(fields[0])
    if not genes:
        die("no gene rows in %s" % path)
    return genes, counts


def read_summary(path):
    names, values = [], {}
    lines = read_lines(path)
    for line in lines[1:]:  # first line is "Status<TAB><BAM path>"
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        names.append(fields[0])
        values[fields[0]] = fields[1]
    if not names:
        die("no status rows in %s" % path)
    return names, values


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=".")
    ap.add_argument("--samples", required=True)
    ap.add_argument("--project-name", required=True)
    ap.add_argument("--out-counts", required=True)
    a = ap.parse_args()

    samples = [s.strip() for s in a.samples.split(",") if s.strip()]
    if not samples:
        die("--samples is empty")
    if len(set(samples)) != len(samples):
        die("duplicate sample ids in --samples: %s" % a.samples)

    first = samples[0]
    genes, _ = read_counts(os.path.join(a.dir, "%s_counts.txt" % first))
    counts_by_sample = {}
    for s in samples:
        g, c = read_counts(os.path.join(a.dir, "%s_counts.txt" % s))
        if s != first:
            if len(g) != len(genes):
                die("sample %s has %d gene rows, first sample %s has %d" % (s, len(g), first, len(genes)))
            for i, (gs, gf) in enumerate(zip(g, genes)):
                if gs != gf:
                    die("sample %s differs from first sample %s at gene row %d: '%s' vs '%s'"
                        % (s, first, i + 1, gs, gf))
        counts_by_sample[s] = c

    out_dir = os.path.dirname(os.path.abspath(a.out_counts))
    with open(a.out_counts, "w", encoding="utf-8") as fh:
        fh.write("\t".join(["Gene"] + samples) + "\n")
        for i, gene in enumerate(genes):
            fh.write("\t".join([gene] + [str(counts_by_sample[s][i]) for s in samples]) + "\n")

    fc_lines = read_lines(os.path.join(a.dir, "%s_featurecounts.txt" % first))
    if len(fc_lines) < 3:
        die("featureCounts table of %s has no data rows" % first)
    with open(os.path.join(out_dir, "gene_length.txt"), "w", encoding="utf-8") as fh:
        for n, line in enumerate(fc_lines[2:], 3):
            fields = line.split("\t")
            if len(fields) < 6:
                die("%s_featurecounts.txt line %d: expected at least 6 columns" % (first, n))
            fh.write("\t".join([fields[0], fields[4], fields[5]]) + "\n")

    status_names = None
    summaries = {}
    for s in samples:
        names, values = read_summary(os.path.join(a.dir, "%s_featurecounts.txt.summary" % s))
        if status_names is None:
            status_names = names
        elif set(names) != set(status_names):
            die("featureCounts summary of sample %s has different status rows than sample %s" % (s, first))
        summaries[s] = values
    with open(os.path.join(out_dir, "featurecounts_summary.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(["Status"] + samples) + "\n")
        for name in status_names:
            fh.write("\t".join([name] + [summaries[s][name] for s in samples]) + "\n")

    print("%s: wrote %s (%d genes, %d samples), gene_length.txt, featurecounts_summary.tsv"
          % (a.project_name, a.out_counts, len(genes), len(samples)))


if __name__ == "__main__":
    main()
