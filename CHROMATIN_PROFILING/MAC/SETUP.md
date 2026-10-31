# SETUP.md — one-time setup of a labmate's Mac for the lab's Docker pipelines

Written for an AI agent (e.g. Claude Code) preparing a labmate's **Mac** to run the
CUT&RUN / chromatin profiling pipeline (`CHROMATIN_PROFILING/MAC`). Do this **once per
computer**. Afterwards the pipeline needs only its own `bash run_cutandrun.sh setup`
(run from the pipeline's folder; it pulls the pipeline's Docker image) and then runs from
the instructions in `README.md` in this folder. Nothing here is specific to one dataset.

The user is usually a wet-lab scientist: explain each step in plain language, ask before
anything that needs an administrator password, a restart, or changes system settings,
and never delete their files.

Origin of this note: the hiccups in the Windows version come from setting up the lab's
first Windows PC (2026-10), and they really happened. The Mac parts are written from
Docker Desktop's documentation and this pipeline's design. They have not yet been
verified on a lab Mac, so expect surprises and report them.

---

## 0. End state (what "set up" means)

1. Docker Desktop installed, its terms accepted, the engine running (`docker info`
   works). Apple-chip Mac: Rosetta 2 installed.
2. Docker's own data (its disk image: images, build cache, and the work volume) on the Mac's
   **internal SSD** with at least 30 GB free, plus room for the work files (about 4.6 times the
   FASTQ total of a dataset, an estimate; see README "Disk space").
3. The pipeline folder cloned (not inside iCloud Drive or a synced Desktop/Documents
   folder), its `setup` run once, and the pipeline's image downloaded
   (`docker image inspect ghcr.io/genesig/omics:<tag>` succeeds; the tag is in
   `docker/image.txt`).
4. A folder for sequencing data, e.g. `~/FASTQ/<run>/`: one subfolder per sequencing
   run; raw FASTQs stay there; projects point at them (they are never moved or copied).
5. A reference folder, e.g. `~/REFERENCES/` (section 4), with the bowtie2 index of the
   combined target + spike-in genome, chromosome sizes, a blacklist, a GTF and a TSS BED.
6. The user knows: keep the Mac plugged in and awake during runs, close memory-hungry
   programs, and how to start and check a run (`README.md`).

---

## 1. Survey first (read-only), then tell the user what you found

**Mac (Terminal):**
```bash
sysctl -n hw.memsize hw.logicalcpu; uname -m; sw_vers; df -h ~ /Volumes/*
diskutil info / | grep -i 'solid state'                  # "Solid State: Yes" = SSD
docker version 2>/dev/null; top -l 1 -o mem | head -20
```
(`uname -m` prints `arm64` on an Apple-chip Mac and `x86_64` on an Intel Mac. Activity
Monitor's Memory tab shows the same as `top`.)

Decide and report:
- **RAM**: 8 GB is the minimum for this pipeline; 16 GB or more is comfortable. With
  Docker Desktop's memory at 6 GB, an 8 GB Mac gives the analysis about 5 GB. A 6 GB Mac
  cannot reach the pipeline's minimum (about 5 GB): the launcher refuses to start a run,
  so tell the user it will not fit.
- **Disks**: which disk is an SSD (the internal disk normally is) and which is a hard
  disk drive (HDD) or a USB drive. Docker's data belongs on the internal SSD (or at least
  on a disk with room). The project folder can be on an HDD or USB drive: the work files
  are on Docker's disk, so only the results are written there (and `intermediate_files/` only with
  `publish_intermediates: true`).
- **Free space**: Docker data several GB (grows with each image pull); the image itself
  is several GB. Each dataset needs about 4.6 times the total size of its FASTQs on Docker's
  disk (the work volume; the launcher computes this before each run; a full dataset can need
  100 GB or more), and room for the results on the project's disk. The raw FASTQs stay where
  they are.
- **Chip**: `arm64` (Apple chip) needs Rosetta 2 (section 2.2); `x86_64` (Intel) does not.
- **Memory hogs**: note big consumers (browsers with many tabs, IDE language servers,
  R or Python sessions with parallel workers). See section 5.

---

## 2. Docker Desktop on the Mac

### 2.1 Install Docker Desktop
The pipeline's `bash run_cutandrun.sh setup` (run from the `MAC/` folder) offers
`brew install --cask docker` when Homebrew exists (it asks first); otherwise it tells the
user to download Docker Desktop for the right chip (Apple chip or Intel) from docker.com.
The user then opens **Docker** from Applications once, approves the helper/password
prompt, accepts the terms (sign-in can be skipped) and waits for "Engine running".
If Docker is installed but not running, `setup` and `run` open it and wait up to 4 minutes.

- **Administrator prompts**: macOS shows them on the user's screen; tell the user to
  approve. A prompt that nobody answers makes the command hang.
- **Restart**: if `docker info` still fails after the first launch, restart the Mac.

### 2.2 Apple-chip Macs: Rosetta 2
The lab images are linux/amd64, which Docker Desktop runs through Rosetta 2 on an Apple
chip. `setup` offers `softwareupdate --install-rosetta --agree-to-license` (it asks first).
Also check Docker Desktop > Settings > General > "Use Rosetta for x86_64/amd64 emulation
on Apple Silicon" (on by default). Without Rosetta the pipeline still runs, much more slowly.

### 2.3 Memory
Docker Desktop > Settings > Resources > **Memory** caps everything Docker's containers may
use. The launcher never changes it. It recommends a value by the Mac's total memory:
below 12 GB, 6 GB; 12 to under 24 GB, 8 GB; 24 GB or more, 16 GB. Set it there, click
*Apply & restart*, and wait for "Engine running".

The analysis needs at least the pipeline's minimum (about 5 GB, `min_memory_gb` in
`docker/requirements.conf`, a placeholder until a pilot run). The budget of a run is the
smallest of 75% of the Mac's memory (at least 3 GB left for macOS), Docker's memory minus
0.5 GB, and 14 GB. So on an 8 GB Mac Docker's memory must be at least 5.5 GB for a run to
start; with 6 GB it gets about 5 GB. `setup` prints what each analysis will get and what
to set if Docker holds it back.

So do not raise Docker's memory "to be safe" by hand: the extra would fill with Docker's
file cache, which macOS gives back slowly, and other programs then slow down. If the user
has changed the setting, leave it unless the analysis is refused or the Mac can spare it.

### 2.4 Files Docker can see
Docker Desktop shares `/Users`, `/Volumes`, `/private`, `/tmp` and `/var/folders` by
default. A FASTQ or reference folder elsewhere needs Docker Desktop > Settings > Resources >
File sharing: the launcher stops with a message until the folder is added. Avoid iCloud-
synced Desktop/Documents folders for big data ("Optimize Mac Storage" can evict files).
External drives under `/Volumes` work (slower if they are USB hard disks).

### 2.5 Double-click launchers
`RUN_ANALYSIS.command` (made by `init`) may be blocked by Gatekeeper the first time:
right-click > Open once. It must stay executable (`chmod +x`).

### 2.6 Notes for the agent's own shells on macOS
- The pipeline's launcher is a bash script that runs on the Mac's own bash (bash 3.2).
  Run it as `bash run_cutandrun.sh ...` if the executable bit is missing.
- The `python3` on macOS may be Apple's stub, which asks to install the Command Line
  Tools. Use the pipeline image's Python instead:
  `docker run --rm --platform linux/amd64 ghcr.io/genesig/omics:<tag> python3 ...`.
- Every `docker run` or `docker pull` you start yourself needs `--platform linux/amd64`
  (the image is amd64; the launcher always passes it).
- Excel or Numbers files: export them as CSV UTF-8 (or convert inside the pipeline image).
- Background jobs that poll a status file can be fooled by the previous run's file for a
  few seconds after `run --detach`; the launcher writes `running` immediately, but prefer
  the launcher's `status` command.

---

## 3. Clone the pipeline and pull its image

- Clone where the user keeps tools, e.g. `~/PIPELINES/`. Not in iCloud Drive or a synced
  Desktop/Documents folder. The labmate gets only this folder with a sparse clone:
  ```bash
  git clone --filter=blob:none --sparse <url> ~/PIPELINES
  cd ~/PIPELINES
  git sparse-checkout set CHROMATIN_PROFILING/MAC
  cd ~/PIPELINES/CHROMATIN_PROFILING/MAC
  ```
  The repository root carries a `.gitattributes` that keeps the Linux scripts LF; it comes
  along with a sparse checkout. Do not edit it.
- On a Mac git keeps the executable bit, so no `core.fileMode` step is needed. If
  `run_cutandrun.sh` or `RUN_ANALYSIS.command` is not executable, run `chmod +x` on it.
- Run the pipeline's setup once, from inside the `MAC/` folder:
  `bash run_cutandrun.sh setup --yes`. It **pulls** the pipeline's image,
  `ghcr.io/genesig/omics:<tag>`, with `docker pull --platform linux/amd64` (the tag is the
  one line of `docker/image.txt`). There is **no image build** on the user's computer: no
  Dockerfile in this folder. It takes several GB of download (keep the Mac online and
  plugged in; corporate proxies and VPNs can break it), then prints the versions of
  Nextflow, R, Python, samtools and bowtie2, and adds `run_cutandrun.sh` to the user's
  PATH (`~/.zshrc`; new Terminal windows only). Run it in the background or with a long
  timeout. `--yes` answers yes to routine questions. `setup` does not change Docker's
  memory setting.
- `run` pulls the image on its own if it is missing, so a failed pull in `setup` is fixed
  by running `setup` again.
- Validation: there is no synthetic test data in this folder yet. Confirm instead that
  `docker image inspect ghcr.io/genesig/omics:<tag>` works and that the `setup` version
  lines print. Full validation of a dataset is the user's first real run (use `status`
  and the QC report, `README.md`).
- The image's tag changes only when `docker/image.txt` changes (`git pull`). A new tag is
  pulled by the next `setup` (or `run`).

---

## 4. Where data and projects go

Example layout on the user's Mac (the layout is a suggestion, not a requirement; the
paths are whatever the user writes in `config.yaml`):

```
~/FASTQ/<run name>/                   raw sequencing runs: one folder per run, the FASTQ files
                                      (R1, R2) and run-level files; read in place, never moved
                                      (an external drive under /Volumes/<name>/ works too)
~/REFERENCES/                         shared reference files, used by several datasets
  bowtie2/t2t_ecoli/                  bowtie2 index prefix of the combined target + spike-in genome
  blacklist/                          blacklist BED
  chr_size/                           chromosome sizes of the combined genome
  genomes/t2t/                        genome files of the t2t build
  annotation/                         gene GTF, TSS BED
~/chromatin/<dataset>/                analysis projects: one folder per dataset (made by `init`)
```

- Raw sequencing runs stay in the FASTQ folder. Several projects can point at the same
  run; the pipeline reads it in place and never moves, copies or renames it. Keep the
  path stable; if the run folder moves, update each project's `config.yaml` and CSVs.
  Network shares written as `smb://...` or `//server/share` are refused by the launcher.
  Mount the share in Finder and use its `/Volumes/<name>/...` path. If Docker cannot read
  that path, copy the run to a local or external drive first.
- Reference files are used in place too. `config.yaml` holds their paths (for example
  `genome_index: /Users/<you>/REFERENCES/bowtie2/t2t_ecoli/<prefix>`).
- Analysis projects: one folder per dataset, created by the pipeline's `init` (e.g.
  `~/chromatin/<dataset>/`), on the internal SSD if possible. If the chosen folder already
  exists and holds anything, list it and ask before using it (`init` adds files and a git
  repository to whatever is there). The project folder must be outside the pipeline folder;
  the launcher refuses otherwise.

---

## 5. Before each long run: free memory and disk

The pipeline caps its own memory, but it can't stop other programs. On the lab's first
Windows PC an IDE's language server grew to 12.5 GB and R sessions held 32 worker
processes (5 GB): with Docker on top, the computer ran out of free memory and stalled.
Before a long run (and if a run or `docker` commands hang), check the top processes
(section 1; Activity Monitor's Memory tab shows the same) and ask the user to close or
restart the big ones. At the start of each run the launcher warns when the Mac's free
memory is below the run's budget.

Before a run:
- check that no other container is running (`docker ps`); the analysis then shares
  Docker's memory with it;
- check the free disk space on Docker's disk (the work volume) and on the project's disk
  against about 4.6 times the FASTQ total on Docker's disk and about 1 time on the project's disk (the launcher prints its numbers at the start of the
  run; the project's disk stops below 20 GB free);
- check Docker's disk is on the internal SSD (or a disk with room); the project folder may be
  on another disk, but not in an iCloud-synced folder;
- confirm with the user that the Mac will stay plugged in with the lid open (the launcher
  keeps it from idle sleep while the analysis runs).

The launcher prints the budget at the start of each run (memory and cores). If the budget
is below the pipeline's minimum, it stops and says why: free memory, close programs, raise
Docker's memory, or use a larger Mac. Do not raise `max_memory_gb` above what the Mac can
spare.

---

## 6. Troubleshooting (symptom → cause → fix)

| Symptom | Likely cause | Fix |
|---|---|---|
| "Docker Desktop cannot see ..." (run stops) | The folder holding an input (or the project) is not shared with Docker | Add the folder under Docker Desktop > Settings > Resources > File sharing, or move the file under `/Users` or `/Volumes`; run again (2.4) |
| Whole Mac stalls during an image pull or a run | Memory pressure: Docker's memory plus other programs | Section 5 (close hogs); Docker's memory per 2.3 |
| `docker` commands hang for minutes; `status` hangs | Docker Desktop starved (low memory, disk saturated) | Free memory; wait; if needed quit Docker Desktop (Docker menu > Quit Docker Desktop, end leftover `com.docker.*` processes in Activity Monitor), open it again from Applications |
| `docker: command not found` in a new Terminal window | Docker's bin folder is not on the PATH yet | Open a new window; the launcher also looks in `/Applications/Docker.app` itself |
| Docker Desktop not starting after `setup` (4 minutes) | First launch and terms not done, approval prompt not answered, or a stuck start | Open Docker from Applications, approve the prompt, accept the terms; quit and reopen it; restart the Mac if it still fails |
| `run` refuses: budget below the pipeline's minimum | Docker's memory too low for the Mac, or `max_memory_gb` set low | Raise Docker's Memory (2.3), free memory, or use a larger Mac; raise `max_memory_gb` only if the Mac can spare it |
| A step ends with exit 137 (killed for memory) | The step needed more than its budget (retried up to three times, never above the budget) | Free memory and run again; set `max_memory_gb` higher only if the Mac can spare it |
| "names a network path" / `smb://` refused | Dataset on a network share written as a network path | Mount the share in Finder and use its `/Volumes/<name>/...` path, or copy the files to a local or external drive |
| Input files reported as not found | A drive is not connected, or the path in `config.yaml`/the CSVs is wrong | Connect the drive, or correct the path |
| `setup` cannot pull the image (network) | No internet, proxy, VPN | Fix connectivity; run `setup` again |
| Analysis much slower than expected | Docker's disk image on a USB hard disk or an HDD, or the project on one with `work_volume: false` | Move Docker's disk image to the internal SSD (Docker Desktop > Settings > Resources > Advanced > Disk image location); with `work_volume: false`, move the project too |
| `RUN_ANALYSIS.command` will not open | Gatekeeper, or the file lost its executable bit | Right-click > Open once (2.5); `chmod +x RUN_ANALYSIS.command` |
| "Permission denied" when running `run_cutandrun.sh` | The executable bit is missing (e.g. after a copy) | `bash run_cutandrun.sh ...`, or `chmod +x run_cutandrun.sh` |

---

## 7. Hand-off to the user

Tell them, briefly: Docker Desktop must be running (it starts with the Mac only if they
enable "Start Docker Desktop when you sign in" in its settings); keep the Mac plugged in
with the lid open during analyses; close big programs before long runs; the launcher never
changes Docker's memory setting, so if it asks for a change they make it in Docker Desktop
(Settings > Resources); raw data lives in the FASTQ folder and is never moved by the
pipeline; how to start a new dataset (`init`, see `README.md`); and that an AI agent can
do all of it for them.
