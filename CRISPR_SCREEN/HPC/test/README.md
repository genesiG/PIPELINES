# Synthetic test project (cluster version)

`make_test_project.py` writes a small, deterministic CRISPR screen (6 samples,
one multiplexed FASTQ of 180k reads, 260 guides incl. non-targeting controls,
6 genes depleted in the treated samples) — enough to run every step of the
pipeline, including demultiplexing, anchor discovery and `mageck test`, in a
few minutes on a laptop. Use it to check a change to the pipeline (or a new
computer's setup) before running real data.

The script needs only Python 3's standard library. Unlike the Docker versions,
nothing runs in a container here: it runs with the `python3` on the login node.
Whether `python3` is installed on the login and compute nodes has not been
checked. If it is missing or too old, load a Python 3 module or activate a conda
environment you have, and use that `python3` instead.

Keep the test project on GPFS (under `~` or `/project/...`), not in `/tmp`:
compute nodes cannot see `/tmp` (see `AGENTS.md`).

**Steps** (from the `HPC` folder, the one that contains `test/`):
```bash
cd /project/rephilllab/CRISPR_SCREEN/HPC   # or the HPC folder of your own clone
out=~/crisprscreen_test
python3 test/make_test_project.py $out
```
This writes `$out/raw/Undetermined_S0_L001_R1_001.fastq.gz` (the multiplexed
FASTQ, outside the project) and `$out/project/` (`samplesheet.csv`,
`library.csv`, `contrasts.csv`, `config.yaml`).

```bash
cd $out/project
run_crisprscreen.sh init
```
`init` works in the current folder and keeps the generator's `config.yaml`
(it only copies the default one when no `config.yaml` exists). The folder does
not need to be empty. `init` also runs `git init` and a first commit here, which
needs `git` on the login node (not checked).

Next, point `fastq_r1` at the generated FASTQ. The generator wrote the
placeholder `RAW_FASTQ_PLACEHOLDER`, which this replaces (GNU `sed`, as on the
cluster):
```bash
sed -i "s#RAW_FASTQ_PLACEHOLDER#$out/raw/Undetermined_S0_L001_R1_001.fastq.gz#" config.yaml
grep '^fastq_r1' config.yaml
```

Then run the pipeline (no arguments) from the project folder:
```bash
run_crisprscreen.sh
```
The launcher checks that the FASTQ exists and submits the Nextflow head job to
LSF (queue `rhel9`). Follow it with:
```bash
bjobs -u $USER
tail -f metadata/nf_logs/nextflow_<timestamp>.out
```
The launcher prints the exact log file names when it submits the job.

**Cheaper checks first** (from `AGENTS.md`, "Validating a change"): these check
the pipeline code, not the full run. Step 2 is a parse and wiring check with
`-preview`. Step 3 is a stub run of the whole pipeline with `-profile local,test
-stub-run`, where every process writes placeholder outputs. Both call Nextflow
directly on the login node with a capped Java heap (`NXF_OPTS`, see `AGENTS.md`),
and are maintainer checks of the pipeline code.

**Reduced resources:** `conf/test.config` lowers the resources of each step for
smoke tests. It is used only when Nextflow is called directly with
`-profile pmacs,test` (or `local,test`). The launcher runs `-profile pmacs`
only, so `run_crisprscreen.sh` in this test uses the normal resources from
`conf/base.config`.

`init` keeps the generator's `config.yaml` (it never overwrites one). Expected:
`Pipeline completed successfully`; every sample demultiplexed with ~30k reads;
an anchor `TTGTGGAAAGGACGAAACACCG` discovered; mapping ~95%; and
`results/mageck_test/treated_vs_baseline/treated_vs_baseline.gene_summary.txt`
with ESS1–ESS6 at the top of the negative selection. Re-running without
changes should report every step as cached.

These expected results were written for the Docker versions. They have not been
checked on the cluster.

**Not verified on the cluster:** `python3` on the login and compute nodes; `git`
for `init`; the generated `config.yaml` as a cluster run (it is not checked
against the cluster launcher); the expected results above; the cached re-run.
The generator's `config.yaml` also has `max_cpus` and `max_memory_gb` keys.
Nothing in `HPC/` reads them, so they have no effect here.
