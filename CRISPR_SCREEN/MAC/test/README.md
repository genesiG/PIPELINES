# Synthetic test project (Docker versions)

`make_test_project.py` writes a small, deterministic CRISPR screen (6 samples,
one multiplexed FASTQ of 180k reads, 260 guides incl. non-targeting controls,
6 genes depleted in the treated samples) — enough to run every step of the
pipeline, including demultiplexing, anchor discovery and `mageck test`, in a
few minutes on a laptop. Use it to check a change to the pipeline (or a new
computer's setup) before running real data.

The script needs only Python's standard library; run it inside the pipeline's
image so the host needs nothing but Docker.

**Mac (Terminal)**, from this `MAC` folder (the one that contains `test/` and
`run_crisprscreen.sh`; start from there, nothing else of the repository is
needed). The Docker image `rephilllab/crisprscreen:0.1.0` must exist first:
`setup` builds it, once per computer (see `AGENTS.md` in this folder):
```bash
bash run_crisprscreen.sh setup --yes
```
Then:
```bash
out=~/crisprscreen_test
docker run --rm -v "$PWD/test:/test:ro" -v "$out:/out" rephilllab/crisprscreen:0.1.0 python3 /test/make_test_project.py /out
bash run_crisprscreen.sh init "$out/project"
sed -i '' "s#RAW_FASTQ_PLACEHOLDER#$out/raw/Undetermined_S0_L001_R1_001.fastq.gz#" "$out/project/config.yaml"
cd "$out/project" && ./utils/run_crisprscreen.sh
```

`init` keeps the generator's `config.yaml` (it never overwrites one). Expected:
`Pipeline completed successfully`; every sample demultiplexed with ~30k reads;
an anchor `TTGTGGAAAGGACGAAACACCG` discovered; mapping ~95%; and
`results/mageck_test/treated_vs_baseline/treated_vs_baseline.gene_summary.txt`
with ESS1–ESS6 at the top of the negative selection. Re-running without
changes should report every step as cached.
