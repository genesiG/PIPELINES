# SETUP.md — one-time setup of a labmate's Mac for the lab's Docker pipelines

Written for an AI agent (e.g. Claude Code) preparing a labmate's **Mac** to run
the lab's containerized pipelines (CRISPR_SCREEN today, RNA-seq, CUT&RUN etc.)
Do this **once per computer**. Afterwards each pipeline only needs its own
`bash run_<pipeline_name>.sh setup` (run from the pipeline's folder; it builds that pipeline's Docker image) and then runs from its runbook (`AGENTS.md` in this
folder). Nothing here is specific to one pipeline.

The user is usually a wet-lab scientist: explain each step in plain language,
ask before anything that needs an administrator password, a restart, or changes
to system settings, and never delete their files.

Origin of this note: the hiccups in it come from setting up the first lab PC (a
Windows 11 machine, 2026-10), and they really happened. The Mac parts are
written from Docker Desktop's documentation and this pipeline's design. They
have not yet been verified on a lab Mac, so expect surprises and report them.

---

## 0. End state (what "set up" means)

1. Docker Desktop installed, its terms accepted, the engine running
   (`docker info` works). Apple-chip Mac: Rosetta 2 installed.
2. Docker's own data (its disk image: images, build cache) on the **internal
   SSD** with ≥ 30 GB free.
3. The pipeline repository cloned (not inside iCloud Drive or a synced
   Desktop/Documents folder), its `setup` run once, and its synthetic test
   project passing.
4. A folder for sequencing data, e.g. `~/FASTQ/<run>/`: one subfolder per
   sequencing run, raw FASTQs stay there; analysis projects point at them (they
   are never moved or copied).
5. The user knows: keep the Mac plugged in and awake during runs, close
   memory-hungry programs, and how to start/check a run (the pipeline runbook).

---

## 1. Survey first (read-only), then tell the user what you found

**Mac (Terminal):**
```bash
sysctl -n hw.memsize hw.logicalcpu; uname -m; sw_vers; df -h ~ /Volumes/*
docker version 2>/dev/null; top -l 1 -o mem | head -20
```
(`uname -m` prints `arm64` on an Apple-chip Mac and `x86_64` on an Intel Mac.)

Decide and report:
- **RAM**: 8 GB or more is comfortable; 6 GB works (the CRISPR pipeline needs
  2.5 GB and runs step by step); less than 6 GB: tell the user it will be slow
  or may not fit.
- **Disks**: the internal disk is normally an SSD (check with
  `diskutil info /` and look for "Solid State: Yes"). External drives may be
  SSDs or hard disk drives (HDD); a USB HDD is the slow case. Docker's data and
  (ideally) the analysis projects belong on the internal SSD. Projects on an
  external HDD work but demultiplexing is several times slower.
- **Free space**: Docker data ~5-20 GB (grows with each pipeline image and
  rebuild); each project needs roughly 1-2x the size of its raw FASTQ (per-sample
  FASTQs + `work/` while a run is going; `work/` is removed after a successful run
  unless `keep_work: true`); the raw FASTQ itself stays in the FASTQ folder.
- **Memory hogs**: note big consumers (browsers with many tabs, IDE language
  servers, R/Python sessions with parallel workers). See section 5.

---

## 2. Docker Desktop on the Mac

### 2.1 Install Docker Desktop
The pipeline's `bash run_crisprscreen.sh setup` (run from the MAC/ folder)
offers `brew install --cask docker` when Homebrew exists; otherwise download
Docker Desktop for the right chip (Apple silicon vs Intel) from docker.com. The
user then opens Docker.app once, approves the privileged-helper/password prompt,
accepts the terms (sign-in can be skipped) and waits for "Engine running".

### 2.2 Apple-chip Macs: Rosetta 2
The lab images are linux/amd64. `setup` offers
`softwareupdate --install-rosetta --agree-to-license`; also enable Docker
Desktop > Settings > General > "Use Rosetta for x86_64/amd64 emulation on Apple
Silicon" when offered. Without Rosetta it still runs, much more slowly.

### 2.3 Memory
Docker Desktop > Settings > Resources > Memory caps everything (default about
half of the RAM). On a 6 GB Mac the default (~3 GB) is enough for the CRISPR
pipeline (it needs 2.5 GB and always leaves 3 GB to macOS). Newer Docker
Desktop versions with the "Docker VMM" virtualization backend reportedly need
≥ 4 GB assigned — if Docker refuses to start, raise it to 4 GB. `setup`
prints what the pipeline would use and what to set if Docker holds it back.

### 2.4 Files Docker can see
Docker Desktop shares `/Users`, `/Volumes`, `/private`, `/tmp` and `/var/folders`
by default; FASTQs elsewhere need Settings > Resources > File sharing. Avoid
iCloud-synced Desktop/Documents for big data ("Optimize Mac Storage" can evict
files). External drives under `/Volumes` work (slower if USB HDD).

### 2.5 Double-click launchers
`RUN_ANALYSIS.command` files may be blocked by Gatekeeper the first time:
right-click > Open once. They must stay executable (`chmod +x`).

---

## 3. Clone the pipeline and build its image

- Clone where the user keeps tools, e.g. `~/PIPELINES/`. Not in iCloud Drive or
  a synced Desktop/Documents folder.
- Only the `MAC/` folder is needed on a labmate's computer, so use a sparse
  clone (it downloads the repository's top-level files and only the folders you
  ask for):
  ```bash
  mkdir -p ~/PIPELINES && cd ~/PIPELINES
  git clone --filter=blob:none --sparse <url>
  cd CRISPR_SCREEN
  git sparse-checkout set MAC
  cd MAC
  ```
- The repository's top-level `.gitattributes` comes with the sparse clone. It
  keeps Linux scripts (`*.sh`, `*.command`, `*.py`, ...) with LF line endings,
  which the container needs. Do not edit it.
- Keep the executable bit. On a Mac git keeps it, so no `core.fileMode` step is
  needed. If `run_crisprscreen.sh` or a `RUN_ANALYSIS.command` is not
  executable, run `chmod +x run_crisprscreen.sh` (and the same for the
  `.command` file).
- Run the pipeline's `setup --yes` once, from inside the MAC/ folder:
  `bash run_crisprscreen.sh setup --yes`. It builds the image (10-30 min,
  ~2-4 GB download, needs internet; corporate proxies/VPNs can break it), prints
  tool versions and adds the command to the user's PATH (new terminals only).
  Run it in the background or with a long timeout. `--yes` answers yes to
  routine questions; it does not change Docker's memory setting (the user raises
  that in Docker Desktop if `setup` asks).
- A transient micromamba crash during the build
  (`Unexpected error 9 on netlink descriptor`) is retried automatically; if the
  whole build fails, simply run `setup` again (finished layers are cached).
- Validate with the pipeline's synthetic test in `test/README.md` (this folder,
  `MAC/test/`): a few minutes; it must end "succeeded" with the expected hits.
  Run it from this MAC folder, as `test/README.md` says (its commands use
  `$PWD/test` and `bash run_crisprscreen.sh`).
- Keep only the pipeline images: helper images pulled for maintenance
  (e.g. BioContainers used to regenerate lock files) can be removed with
  `docker rmi`.

---

## 4. Where data and projects go

- Raw sequencing runs: `~/FASTQ/<run name>/` or an internal-SSD folder of the
  user's choice — one subfolder per run, with the multiplexed
  `Undetermined*_R1_001.fastq.gz` (and R2), plus run-level files (full sample
  sheet, run QC). Several labmates' projects can point to the same run; the
  pipelines read it in place and never move, copy or rename it. Keep the path
  stable; if the run folder moves, update each project's `config.yaml`. Network
  shares can't be read by Docker Desktop: copy the run to a local or external
  drive first.
- Analysis projects: one folder per screen, created by the pipeline's `init`
  (e.g. `~/CRISPR_screens/<name>/`), on the internal SSD if possible. If the
  chosen folder already exists and holds anything, list it and ask before using
  it (`init` adds files and a git repository to whatever is there).

---

## 5. Before each long run: free memory

The pipelines cap their own memory, but they can't stop other programs. On a
lab machine an IDE's language server grew to 12.5 GB and R sessions held 32
worker processes (5 GB): with Docker on top, the computer ran out of free
memory and stalled. Before a long run (and if a run or `docker` commands hang),
check the top processes (section 1; Activity Monitor's Memory tab shows the same)
and ask the user to close or restart the big ones. The CRISPR launcher warns
when free memory is below the run's budget.

---

## 6. Troubleshooting (symptom → cause → fix)

"Scope" says whether the row is specific to the Mac or generic (applies to any
computer running these pipelines, with the Mac wording used here).

| Symptom | Likely cause | Fix | Scope |
|---|---|---|---|
| "mounts denied" / path not shared | FASTQ outside Docker's shared folders | 2.4 | Mac |
| Very slow on Apple chip | Rosetta missing | 2.2 | Mac |
| `docker` commands hang for minutes; `docker info` hangs | Docker Desktop starved (low memory, disk saturated) | Free memory (section 5); wait; if needed quit Docker Desktop and start it again | Generic |
| Docker Desktop starved or not starting | Low memory, or Docker Desktop stuck from earlier | Quit Docker Desktop (Docker menu > Quit Docker Desktop; end any leftover `com.docker.*` processes in Activity Monitor), open it again from Applications | Mac |
| `docker info` fails after install | First launch/terms not done, or restart pending | Open Docker.app from Applications, approve the helper prompt, accept terms; restart the Mac if it still fails | Generic (Mac wording) |
| Image build fails (network) | No internet, proxy, VPN | Fix connectivity; rerun `setup` | Generic |
| Analysis much slower than expected | Project or FASTQ on an external USB HDD | Move the project/data to the internal SSD or a fast external SSD | Mac |

---

## 7. Hand-off to the user

Tell them, briefly: Docker Desktop must be running (it starts with the computer
only if they enable "Start Docker Desktop when you sign in" in its settings);
keep the Mac plugged in and awake during analyses; close big programs before
long runs; raw data lives in the FASTQ folder and is never moved by the
pipelines; how to start a new screen (the pipeline runbook's `init`), and that
an AI agent can do all of it for them.
