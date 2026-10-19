# AGENTS.md — CRISPR_SCREEN pipeline, cluster (HPC) version (Rephill Lab)

> **Running an analysis on a personal computer (Windows / Mac)?** Use the
> Docker versions: the `WINDOWS/` or `MAC/` folder of the full repository,
> each with its own runbook (`AGENTS.md`). Nothing in this file is needed for
> that. A single-folder clone (sparse checkout of `HPC/`) does not contain
> them.
>
> **Repository layout:** `HPC/` (this folder) is the cluster (PMACS/LSF)
> version described below. In the full repository, `WINDOWS/` and `MAC/` are
> self-contained Docker ports of it: same DAG, modules and `bin/` scripts, one
> Docker image instead of per-step Apptainer containers, local executor instead
> of LSF. A change to the pipeline's logic (`main.nf`, `modules/`, `bin/`)
> usually belongs in all three. `tools/sync_os_folders.sh` (full repository)
> keeps WINDOWS/ and MAC/ identical and checks the `bin/` and agent kits. The
> porting rules for the local Docker versions are in section B of
> `WINDOWS/AGENTS.md` / `MAC/AGENTS.md` (full repository).
>
> **First-time setup on the cluster:** `SETUP.md`, once. **Subagents:**
> `.claude/agents/` in this folder is the same kit as `agent/.claude/agents/`
> (which `init` copies into every project); who does what is in `CLAUDE.md`.

This is the **shared, lab-wide** Nextflow pipeline for pooled CRISPR screens on
UPenn PMACS (LSF). The repository is cloned once to
`/project/rephilllab/CRISPR_SCREEN`; every lab member's projects run the code
in this `HPC/` folder, through `/project/rephilllab/bin/run_crisprscreen.sh` →
`/project/rephilllab/CRISPR_SCREEN/HPC/run_pipeline.sh` (a symlink; after
moving the cluster version into `HPC/`, re-point it once with
`ln -sfn /project/rephilllab/CRISPR_SCREEN/HPC/run_pipeline.sh /project/rephilllab/bin/run_crisprscreen.sh`).

This was the lab's first Nextflow port and is the reference implementation
for porting others (RNA-seq, CUT&RUN, ...). See
`/project/rephilllab/NEXTFLOW_PORTING_GUIDE.md` for the conventions,
PMACS/Nextflow gotchas, and validation methodology distilled from building
this one — read it before starting a new port.

## Before you change anything

- **Only modify this repository when the user explicitly asks to change the
  lab-wide pipeline.** If they want to change *their own* analysis or plots,
  that belongs in their project directory as a `config.yaml` override (see
  `agent/AGENTS.md` and the README's "Customizing" section), not here.
- If you're helping with a **project** (a directory with `config.yaml`,
  `samplesheet.csv`, `results/`), follow that project's `AGENTS.md` instead.
- Changes here affect everyone's next run. Changing a module's inputs or
  script invalidates that step's cache for every project, so it re-runs once.
  Say so when you report a change.

## Layout

```
run_pipeline.sh        init (scaffold a project) / run (check FASTQ paths, bsub the Nextflow head job)
main.nf                the DAG; stepScripts() resolves each step's script (bin/ or the project's replacement)
nextflow.config        params (defaults), profiles, provenance reports -> metadata/pipeline_info/
conf/base.config       per-process resources, retry policy, shell = bash -euo pipefail
conf/pmacs.config      LSF executor + Apptainer (incl. the FIPS bind mount) + shared container cache
conf/test.config       small resources for smoke tests;  conf/conda.config: unused fallback
modules/local/*.nf     one process per step; functions.nf: runner() / extraArgs() shared by all of them
bin/                   every step's script (.sh wrappers around the tools, Python/R helpers) + my_theme.R
examples/              files `init` copies into a project (config.yaml template, *_example.csv, project.gitignore)
agent/                 AGENTS.md / CLAUDE.md, and .claude/ (agent definitions), that `init` copies into every project; init also copies bin/plot_barcode_position.R to utils/ (if absent)
.claude/agents/        subagent definitions for working on this folder (same kit as agent/.claude/agents/; table in CLAUDE.md)
SETUP.md               first-time setup on the cluster (once per user; the lab's shared steps once)
test/                  make_test_project.py + README.md: a small synthetic screen for validating a change
```

## Conventions (keep them when adding or changing steps)

- **Every step's work is a script in `bin/`, passed in as a task input** (not
  called off PATH, and no tool commands inline in a module). A project
  replaces a script with a `config.yaml` key named after it (file name without
  extension); `stepScripts()` in `main.nf` builds that mapping from whatever is
  in `bin/`, so a new script is replaceable with no param declaration. Modules
  call it as `${runner(script)} "${script}" ... ${extraArgs('<name>')}` (from
  `modules/local/functions.nf`). A script's header must document the
  arguments it receives and the files it must write, since that's the contract
  a replacement has to meet. `.sh` scripts pass unrecognized arguments on to
  their tool. A new script needs a row in the README's customization table
  and in `agent/AGENTS.md`. Its name must not collide with an existing param
  (e.g. don't name one `library.py`).
- **Compute and plotting are separate processes.** A plotting step reads only
  tables that an upstream step wrote, so plot changes never recompute
  anything. Examples: `GUIDE_REPRESENTATION_WRANGLE` → `GUIDE_REPRESENTATION_PLOT`,
  and `DEMULTIPLEX_*` (writes the barcode-position CSV) → `PLOT_BARCODE_POSITION`.
  `GUIDE_REPRESENTATION_REPORT` extends this one step further: it only reads
  WRANGLE's tables and PLOT's SVGs (joined by library-name prefix, same
  `failOnMismatch`/`failOnDuplicate` rule as everywhere else) and writes one
  HTML page — it never recomputes a table or re-renders a plot itself.
  `GUIDE_REPRESENTATION_COMPARE_REPORT` (only with >1 library) is the one
  exception to the split: it builds its own comparison plots *and* the report
  in one script/process. It's still downstream-only (reads every library's
  WRANGLE tables, `.collect()`ed, never a count table directly), so it's still
  cheap to rebuild — the split exists to avoid re-running expensive
  computation on a plot-only change, and there's no expensive computation
  here to separate from the plotting.
- **Plot scripts** take `--key=value` args and `--theme=<file>` (defining
  `my_theme`). They fall back to `my_theme.R` next to themselves via
  Rscript's `--file=` argument. Write SVG with base `svg()`; `ggsave(device =
  "svg")` needs svglite, which `rocker/tidyverse:4.3` lacks.
- **Every process publishes its own exact command as `command.sh`**, next to
  its `.log` file, so a labmate or agent can see precisely what ran without
  digging into Nextflow's ephemeral `work/` dirs — this is what the legacy
  pipeline's own per-step batch files did. Mechanism: `conf/base.config` sets
  a global `afterScript = 'cp .command.sh command.sh ...'` (runs after the
  real script *or* the stub, regardless of exit status — the same line
  covers every process with no per-module Groovy needed for the copy
  itself). Each process then just needs to (1) add `path "command.sh", emit:
  command_sh` to its `output:`, and (2) add its own `publishDir ...,
  pattern: "command.sh", saveAs: { "<meaningful-name>_command.sh" }`. **The
  existing publishDir line(s) almost always need a `pattern:` added too** —
  most of them had none (default: publish every declared output), which
  would otherwise also sweep up the generic `command.sh` unrenamed and
  collide across tasks; see any module for the exact pattern. Only add a
  bare `pattern: "..."` (not a closure) when it uses `params.*` — a
  **task-input variable in `pattern:`/the publishDir path must be inside a
  closure** (`pattern: { "${prefix}..." }`, `publishDir { "...${name}" }`),
  the same rule as `saveAs:` — a real bug this session used a bare string.
  Publishing only happens for a *successful* task (standard Nextflow
  behavior, not fixable) — a failed task's `command.sh` still exists (
  `afterScript` ran) but only in its work dir, which is why step 2 of "when
  a run fails" still points there first.
- **Small-input steps downstream of deterministic tools use `cache 'deep'`**
  (content-hashed inputs), so an identical re-generated upstream output
  doesn't force a re-run. Don't use it on steps whose inputs are FASTQs:
  hashing multi-GB files is slow.
- **Every process has a `stub:` block** that writes placeholder outputs
  matching its `output:` declarations. Keep them when editing a process, and
  add one to any new process: `-stub-run` (validation step 3 below) depends on
  them. If `main.nf` recovers a channel's key from an output *filename*
  (e.g. `guide_tables_by_pool`'s `anchor.name - '_counts_long.csv'`), the stub
  must write that same real filename, not a generic placeholder — a stub
  named `${prefix}_stub.png` instead of `${prefix}_ecdf.svg` passed `-stub-run`
  on its own process but broke the *next* process's join with a null-pointer
  error, since nothing existed to recover the key from (a real bug found this
  way — see `GUIDE_REPRESENTATION_PLOT`'s stub).
- **Channels keyed by library name** (`join(..., failOnMismatch: true,
  failOnDuplicate: true)`): never pair per-library outputs by channel order.
  Without `failOnMismatch`, a key mismatch silently drops that library and the
  run still "succeeds".
- **`splitCsv(header: true, quote: '"')`** for any CSV with quoted fields.
  `validated_contrasts.csv` has quoted comma lists, and without `quote` the
  fields shift.
- **Exit code 3 = known data problem.** `errorStrategy` finishes on it without
  retrying. Use it with a plain-language message saying what the user should
  fix, as `MAGECK_COUNT` does for samples with zero mapped reads.
- **Directory convention for projects:** input CSVs + `config.yaml` at the
  project root; `original_data/fastq/` (raw FASTQs stay where config.yaml points); `metadata/`
  (pipeline-generated run bookkeeping); `results/<step>/`, with tool-format
  conversions in `results/intermediate_files/` (default `intermediate_dir`).
  Output filenames are prefixed with `params.project_name`.
- **`METHODS.md` must stay true.** `bin/write_methods.py` (process `METHODS`,
  last in the DAG) writes the run's Materials and Methods from `runInfo()`
  in `main.nf` plus the steps' own evidence: tool logs (cutadapt and MAGeCK
  record their exact command line and version), anchors, libraries,
  `*.sample_reads.tsv`, contrasts, and `SOFTWARE_VERSIONS`' probe of every
  container image. Tool flags are read from the logs, so changing a flag in a
  `.sh` script is described correctly, but the prose explaining *what a step
  does* (algorithms, QC metric definitions, default thresholds) lives in
  `write_methods.py`. When you change a step's behaviour, defaults, outputs
  or file names, or add a step or an R package, update its section there
  (and `R_PACKAGES` in `bin/software_versions.sh`), then check the output on
  real evidence, not only with `-stub-run` (whose METHODS.md is placeholder
  text). A new step's evidence goes into the `evidence` channel in `main.nf`;
  its file names must be unique across all evidence (one staging directory).
- Keep `README.md`, `examples/config.yaml` and `agent/AGENTS.md` consistent
  with any change to params, outputs or commands.

## Known traps (each of these broke a real run)

- **Embedding several R `svg()` outputs inline in one HTML page.** Cairo (the
  device backend) numbers each file's own `<symbol>`/`<clipPath>` ids from
  scratch (`glyph0-0`, `clip1`, ...). Two independently-rendered SVGs in the
  same page therefore reuse the same ids for different glyph shapes, and the
  browser resolves every `<use xlink:href="#glyph0-0">` to whichever plot's
  definition wins — text renders with other plots' letters ("scrambled"
  titles/labels; a real bug, found by the user, not caught by any of the
  validation steps below since the HTML was well-formed and parsed fine).
  Both `guide_representation_report.R` and `guide_representation_compare_report.R`
  fix this with a `namespace_svg_ids()` helper that prefixes every id (and
  every `href="#id"`/`url(#id)` reference to it) with something unique to that
  plot before embedding it. Any future script that embeds more than one
  R-rendered SVG in one page needs the same fix.
- **FIPS:** every Apptainer exec on PMACS compute nodes panics unless
  `/usr/lib64/ossl-modules` is bind-mounted (`conf/pmacs.config`
  `runOptions`). Don't remove it.
- **Nextflow 26 strict syntax:** no bare statements at top level of `main.nf`
  or config (put logic in functions or `workflow {}`); `publishDir` paths that
  use task inputs must be closures; `workflow.onComplete` lives in
  `nextflow.config` (and uses `println`, since `log` is unbound there). A
  closure stored in a local variable can't be called like a function
  (`def f = { ... }; f(x)` fails to parse: "`f` is not defined") — write a
  top-level `def` function instead.
- **`publishDir`'s `pattern:` does NOT support a closure at all** — unlike the
  destination path and `saveAs:`, which both accept one. `pattern: {
  "${prefix}_qc_report.html" }` silently matches nothing: no error, no
  warning, the task still "succeeds", the file just never appears under
  `results/`. A real bug this session shipped this way in an already-committed
  batch-file-logging change, in both `guide_representation_report.nf` (the
  main QC report HTML) and `mageck_test.nf` (the actual `mageck test`
  outputs) — only `command.sh` published, the real output silently didn't.
  Not caught by that change's own regression check because it only grepped
  for `command.sh`'s presence, never re-verified the main output still
  published. Fix: never put a closure in `pattern:`; instead give the process
  two `publishDir` lines with no `pattern:` restriction at all, one for
  `results/<dir>` with `saveAs: { fn -> fn == 'command.sh' ? null : fn }`
  (publishes everything except `command.sh`) and one for `results/<dir>/log`
  with a bare-string `pattern: "command.sh"` and its own `saveAs:` to rename
  it — see `mageck_test.nf` or `guide_representation_report.nf`. After
  touching any `publishDir`/`pattern:` line, re-check with `-stub-run` (or a
  real run) that the process's ACTUAL declared outputs appear under
  `results/`, not just that `command.sh` does.
- **Pipelines through `tee`:** `pipefail` is set globally. Without it, a
  crashing tool piped into `tee` exits 0 and fails later as a misleading
  "missing output file".
- **`cpus` + `clusterOptions` together** produce a duplicate
  `-R "span[hosts=1]"` that LSF rejects. `DEMULTIPLEX_*` therefore sets only
  `clusterOptions`, and runs `cutadapt --cores 0`.
- **`Sys.which()` returns `""` for non-executable files.** Don't use it to
  locate helper files. New `bin/` scripts still need `chmod +x`.
- **MAGeCK normalizes by each sample's total mapped reads**, even with
  `--norm-method median`, so a zero-mapping sample crashes `mageck count`. It's
  handled with an exit-3 message; don't try to "fix" it by altering counts.
- **Every sample against every library (no samplesheet `library` column,
  >1 library) with anchor discovery:** each library's `GUIDE_ANCHOR_TRIM`
  writes the same `<sample>_trimmed.fastq.gz` names, so `MAGECK_COUNT` must
  get only its own library's trimmed files (they travel in its input tuple,
  joined by library name). Pooling every library's trimmed FASTQs into one
  input list made `MAGECK_COUNT` fail with "input file name collision" — a
  real bug, present until the `-stub-run` of this mode caught it. The
  published copies go to `anchor_trimmed_fastq/<lib>/` for the same reason.
- **Anchor discovery must be checked against the library.** Consensus
  extension alone absorbed the guide's first base into the anchor for a
  library whose 20-mers all start with the vector `G` (`G` + 19 nt; real
  case, AK_SCREEN_2026): the anchor came out 22 nt instead of 21, every
  trimmed read was shifted by one, and nothing mapped. `discover_anchor.py
  --library` now trims candidates back from the 3' end and keeps the one
  followed by an exact library guide in the most reads.
- **`mageck count` has no mismatch-tolerance option** — verified directly
  against `mageck count --help` (0.5.9.5), not assumed. Its guide matching is
  always an exact hash lookup on the trimmed read; `--sgrna-len` only tells
  it how to format `--unmapped-to-file`'s output (the library file itself
  determines the actual matching length/content). The former
  `target_mismatches` key (removed) only named output files and never reached
  mageck — don't assume it changed matching behavior. mageck
  count's own mapping rate already *is* "% exact matches"; the near-miss/
  far-miss classification (`guide_representation_wrangle.R`) is built
  separately, from `--unmapped-to-file`'s output.
- The shared scratch/work dir is `/project/rephilllab/nf_work/<user>/<project>`.
  The head job deletes the work folder and Nextflow's `.nextflow/` cache after a
  successful run unless `config.yaml` has `keep_work: true` (default `false`); failed
  or stopped runs keep both, so a re-run resumes. Re-running a finished project recomputes
  everything (~10-20 min for a full lane). `publish_mode` is `copy`, so nothing
  in `results/` points into the work dir and the cleanup is safe.
  Session-scoped `/tmp` paths are not visible from compute nodes, so test
  inputs must live on GPFS (`/home`, `/project`).

## Validating a change (cheapest first; don't re-run labmates' data)

`test/README.md` in this folder describes the synthetic test project
(`test/make_test_project.py`): a small screen that runs every step in a few
minutes. Use it as the input for steps 3 and 4 below, and for step 5 when the
user asks for a full run.

1. `bash -n run_pipeline.sh` for shell edits.
2. **Parse + DAG wiring on the login node** (no tasks run). Java needs a
   capped heap there:
   ```bash
   export NXF_OPTS="-Xmx384m -Xms64m -XX:MaxMetaspaceSize=160m"
   export NXF_HOME=/project/rephilllab/nf_work/.nextflow_home
   # in a scratch dir with dummy samplesheet.csv / library.csv / reads.fastq.gz / config.yaml:
   /project/rephilllab/CRISPR_SCREEN/HPC/.nextflow-bin/nextflow -log p.log run /project/rephilllab/CRISPR_SCREEN/HPC \
       -profile pmacs -params-file config.yaml -preview
   ```
   `-preview` doesn't execute script blocks, so it won't catch runtime errors
   inside them.
3. **Stub run of the whole DAG on the login node**: every process runs its
   `stub:` block (placeholder files, the local executor, no containers, no
   LSF), so channel wiring, joins, and startup validation execute for real, in
   seconds. Use a scratch project with dummy inputs (per-sample FASTQs can be
   one-line placeholders in `fastq_dir` mode):
   ```bash
   .../nextflow -log nf.log run /project/rephilllab/CRISPR_SCREEN/HPC -profile local,test \
       -params-file config.yaml -stub-run -ansi-log false
   ```
   Then check `results/` (e.g. count-table headers, `validated_contrasts.csv`).
   Test the error paths too; each should fail at startup, or with exit 3, with
   a clear message.
4. **Run the changed script alone in its container**, submitted with `bsub`
   (small resources), on real or synthetic inputs, with `bin/` bound the way
   Nextflow mounts it:
   ```bash
   apptainer exec -B /usr/lib64/ossl-modules -B $(pwd) -B /project/rephilllab/CRISPR_SCREEN/HPC/bin:/nf_bin \
       /project/rephilllab/nf_containers/<image>.img Rscript /nf_bin/<script>.R --...
   ```
   Compare against previous outputs (for example, `diff` the tables) when the
   change should preserve results.
5. **A full pipeline run only when the user asks for it.** Use a throwaway
   project on GPFS with synthetic or copied data, never a labmate's project.
   Reference result: for the `AV_CRISPRlib` sample, `mageck count` must give
   `29d909723348548d420e6581ac937fd1` (md5 of its `.count.txt`), identical
   to the pre-Nextflow pipeline.

Never run Nextflow tasks, analysis scripts, or anything memory-heavy directly
on the login node (~2 GB cap).

## Git

The repository is git (its root is the `CRISPR_SCREEN` folder, one level above
this one). Commit each coherent change with a message that says what changed,
why, and how it was validated. Check `git status` before and after.

- Shell scripts (`*.sh`, including `run_pipeline.sh` and `bin/*.sh`) must stay
  executable: `git add --chmod=+x <file>` for a new one. They must stay LF
  (`.gitattributes`); a CRLF script fails inside Apptainer.
- `.bat` and `.ps1` files (the Windows launchers) keep CRLF; they live in
  `WINDOWS/`, not here.
- Don't commit `.nextflow-bin/`, `.nextflow*`, work dirs, project folders, or
  data (results, FASTQs, `*.sif`).
