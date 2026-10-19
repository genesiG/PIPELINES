#!/usr/bin/env python3
"""make_test_project.py — a small synthetic CRISPR screen, for testing the pipeline.

Writes, under OUTDIR:
  raw/Undetermined_S0_L001_R1_001.fastq.gz   one multiplexed single-end FASTQ
                                              (outside the project, so staging is exercised)
  project/samplesheet.csv, library.csv, contrasts.csv, config.yaml

Read structure (what the real screens look like):
  [0-3 nt stagger][8 nt sample barcode][TTGTGGAAAGGACGAAACACCG anchor][20 nt guide][vector tail]
6 samples (baseline / treatment, 3 replicates each), 60 genes x 4 guides + 20
non-targeting controls; genes ESS1..ESS6 are depleted in treatment, so
`mageck test` has real hits to find. ~1-2% of reads carry a sequencing error
in the guide (near-miss) and ~2% are junk (far-miss). Deterministic (seeded).

Stdlib only. Run it inside the pipeline's image (no Python needed on the host):
  docker run --rm -v "<repo>/test:/test:ro" -v "<OUTDIR>:/out" rephilllab/crisprscreen:0.1.0 \
      python3 /test/make_test_project.py /out
"""
import csv
import gzip
import os
import random
import sys

ANCHOR = "TTGTGGAAAGGACGAAACACCG"
TAIL = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC"
READ_LEN = 75
SAMPLES = [  # barcode_id, barcode, sample_name
    ("BC1", "GCCATCTT", "baseline_rep1"),
    ("BC2", "GCTGGATT", "baseline_rep2"),
    ("BC3", "TAACTCGG", "baseline_rep3"),
    ("BC4", "ATGCAGTC", "treated_rep1"),
    ("BC5", "CGTATGCA", "treated_rep2"),
    ("BC6", "TCGACTAG", "treated_rep3"),
]
READS_PER_SAMPLE = 30000


def main(outdir: str) -> None:
    rng = random.Random(20261015)
    project = os.path.join(outdir, "project")
    raw = os.path.join(outdir, "raw")
    os.makedirs(project, exist_ok=True)
    os.makedirs(raw, exist_ok=True)

    def rand_seq(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    # Library: 60 genes x 4 guides + 20 non-targeting controls, unique 20-mers.
    seen, library = set(), []
    genes = [f"ESS{i}" for i in range(1, 7)] + [f"GENE{i}" for i in range(1, 55)]
    for g in genes + ["NTC"] * 5:
        for k in range(4):
            s = rand_seq(20)
            while s in seen:
                s = rand_seq(20)
            seen.add(s)
            gid = f"{g}.{k + 1}" if g != "NTC" else f"NTC.{len([x for x in library if x[0] == 'NTC']) + 1}"
            library.append((g, gid, s))
    with open(os.path.join(project, "library.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "gRNA_ID", "sequence"])
        w.writerows(library)

    with open(os.path.join(project, "samplesheet.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode_id", "barcode", "sample_name"])
        w.writerows(SAMPLES)

    with open(os.path.join(project, "contrasts.csv"), "w", newline="") as fh:
        fh.write("treatment,control\ntreated,baseline\n")

    # Guide abundance: log-normal-ish skew; ESS* depleted 8x in treated samples.
    base = [rng.lognormvariate(0, 0.5) for _ in library]
    records = []
    for bc_id, bc, name in SAMPLES:
        treated = name.startswith("treated")
        w = [b / (8.0 if treated and g.startswith("ESS") else 1.0) for b, (g, _, _) in zip(base, library)]
        cw, total = [], 0.0
        for x in w:
            total += x
            cw.append(total)
        for i in range(READS_PER_SAMPLE):
            r = rng.random()
            if r < 0.02:                                   # junk read: far-miss
                insert = ANCHOR + rand_seq(20)
            else:
                guide = rng.choices(library, cum_weights=cw)[0][2]
                if r < 0.035:                              # one sequencing error in the guide: near-miss
                    p = rng.randrange(20)
                    guide = guide[:p] + rng.choice([c for c in "ACGT" if c != guide[p]]) + guide[p + 1:]
                insert = ANCHOR + guide
            seq = (rand_seq(rng.randrange(4)) + bc + insert + TAIL + rand_seq(30))[:READ_LEN]
            records.append((f"@SYN:{name}:{i}", seq))
    rng.shuffle(records)
    fq = os.path.join(raw, "Undetermined_S0_L001_R1_001.fastq.gz")
    with gzip.open(fq, "wt") as fh:
        for rid, seq in records:
            fh.write(f"{rid}\n{seq}\n+\n{'I' * len(seq)}\n")

    with open(os.path.join(project, "config.yaml"), "w") as fh:
        fh.write(
            "# Synthetic test project (test/make_test_project.py)\n"
            "project_name: SYNTHETIC_TEST\n"
            "samplesheet: samplesheet.csv\n"
            "library: library.csv\n"
            "# Set by the test: where the generator wrote the FASTQ, as the HOST sees it.\n"
            "fastq_r1: RAW_FASTQ_PLACEHOLDER\n"
            "fastq_r2: null\n"
            "fastq_dir: null\n"
            "paired_end: false\n"
            "sgrna_len: 20\n"
            "contrasts: contrasts.csv\n"
            "max_cpus: null\n"
            "max_memory_gb: null\n"
            "original_data_dir: ./original_data\n"
            "intermediate_dir: ./results/intermediate_files\n"
            "outdir: ./results\n"
        )
    print(f"Wrote {len(records)} reads to {fq}")
    print(f"Project files in {project}: samplesheet.csv library.csv contrasts.csv config.yaml")
    print("Set fastq_r1 in config.yaml to the FASTQ's path on the host before running.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
