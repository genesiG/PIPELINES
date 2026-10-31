# Synthetic test project (Docker versions)

`make_test_project.py` writes a small, deterministic CRISPR screen (6 samples,
one multiplexed FASTQ of 180k reads, 260 guides incl. non-targeting controls,
6 genes depleted in the treated samples) — enough to run every step of the
pipeline, including demultiplexing, anchor discovery and `mageck test`, in a
few minutes on a laptop. Use it to check a change to the pipeline (or a new
computer's setup) before running real data.

The script needs only Python's standard library; run it inside the pipeline's
image so the host needs nothing but Docker.

**Windows (PowerShell)**, from the `WINDOWS` folder (the one that contains `test\`).
The Docker image `rephilllab/crisprscreen:0.1.0` must exist first: `setup` builds
it, once per computer (10–30 minutes the first time, see `AGENTS.md`):
```powershell
.\run_crisprscreen.bat setup --yes
```
Then:
```powershell
cd "<path to this WINDOWS folder>"
$out = "$env:USERPROFILE\crisprscreen_test"
docker run --rm -v "${PWD}\test:/test:ro" -v "${out}:/out" rephilllab/crisprscreen:0.1.0 python3 /test/make_test_project.py /out
.\run_crisprscreen.bat init "$out\project"
# config.yaml: point fastq_r1 at the generated FASTQ (outside the project, so reading a FASTQ from outside the project is tested too)
(Get-Content "$out\project\config.yaml") -replace 'RAW_FASTQ_PLACEHOLDER', "$out\raw\Undetermined_S0_L001_R1_001.fastq.gz" | Set-Content "$out\project\config.yaml"
cd "$out\project"; .\utils\run_crisprscreen.bat run --detach
.\utils\run_crisprscreen.bat status
```

`init` keeps the generator's `config.yaml` (it never overwrites one). Expected:
`Pipeline completed successfully`; every sample demultiplexed with ~30k reads;
an anchor `TTGTGGAAAGGACGAAACACCG` discovered; mapping ~95%; and
`results/mageck_test/treated_vs_baseline/treated_vs_baseline.gene_summary.txt`
with ESS1–ESS6 at the top of the negative selection. Re-running without
changes should report every step as cached.
