# SETUP.md — one-time setup of a labmate's Windows PC for the lab's Docker pipelines

Written for an AI agent (e.g. Claude Code) preparing a labmate's **Windows PC** to run
the RNA-seq pipeline (`RNASEQ\WINDOWS`). Do this **once per computer**. Afterwards the
pipeline needs only its own `run_rnaseq setup` (which pulls the pipeline's Docker image)
and then runs from the instructions in `README.md` in this folder. Nothing here is specific
to one dataset.

The user is usually a wet-lab scientist: explain each step in plain language, ask before
anything that needs an administrator prompt, a restart, or changes system settings, and
never delete their files.

The Docker steps below are the same as for the lab's other Docker pipelines (the first of
them was set up on the lab's first Windows PC, 2026-10). The RNA-seq specific parts are the
inputs and the reference files.

---

## 0. End state (what "set up" means)

1. Docker Desktop installed, its terms accepted, the engine running (`docker info`
   works), on WSL 2.
2. Docker's own data (its disk image: images, build cache, and the work volume) on an
   **internal SSD** with at least 30 GB free, plus room for the work files (see README
   "Disk space"; the RNA-seq figure is not measured yet).
3. The page file (virtual memory) on an SSD, not on a hard disk drive.
4. The pipeline folder cloned (not inside OneDrive-synced folders), its `setup` run once,
   and the pipeline's image downloaded (`docker image inspect ghcr.io/genesig/omics:<tag>`
   succeeds; the tag is in `docker\image.txt`).
5. A folder for sequencing data, e.g. `D:\FASTQ\<run>\`: one subfolder per sequencing
   run; raw FASTQs stay there; projects point at them (they are never moved or copied).
6. A reference folder, e.g. `D:\REFERENCES\` (section 4), with the HISAT2 index of the
   genome (prefix) and the gene annotation GTF.
7. The user knows: keep the computer plugged in and awake during runs, close
   memory-hungry programs, and how to start and check a run (`README.md`).

---

## 1. Survey first (read-only), then tell the user what you found

**Windows (PowerShell):**
```powershell
$cs = Get-CimInstance Win32_ComputerSystem; $os = Get-CimInstance Win32_OperatingSystem
"RAM {0:N1} GB, free {1:N1} GB, logical cores {2}" -f ($cs.TotalPhysicalMemory/1GB), ($os.FreePhysicalMemory/1MB), $cs.NumberOfLogicalProcessors
$os.Caption + ' ' + $os.Version
Get-PhysicalDisk | Select FriendlyName, MediaType, BusType, @{n='GB';e={[int]($_.Size/1GB)}}
Get-Partition | ? DriveLetter | Select DriveLetter, DiskNumber        # which drive letter is which disk
Get-PSDrive -PSProvider FileSystem | Select Name, @{n='FreeGB';e={[int]($_.Free/1GB)}}
Get-CimInstance Win32_PageFileUsage | Select Name, AllocatedBaseSize   # where virtual memory lives
(Get-CimInstance Win32_Processor).VirtualizationFirmwareEnabled        # must be True for WSL 2
wsl --status; docker version                                           # may fail: not installed yet
Get-Process | Sort WorkingSet64 -Desc | Select -First 10 Name, @{n='MB';e={[int]($_.WorkingSet64/1MB)}}
```

Decide and report:
- **RAM**: the pipeline's minimum budget and the most it can use are in `docker\requirements.conf` (`min_memory_gb`,
  `max_useful_memory_gb`). Compare them with the budget the computer gets (75% of RAM, leaving 3 GB to Windows, and Docker's
  own limit). Below the minimum the launcher refuses to start a run; tell the user if the computer will not fit the full data.
- **Disks**: which drive is an SSD (`MediaType SSD`/NVMe) and which is a hard disk drive
  (`HDD`). Docker's data (its disk image, which holds the work volume) and the page file
  belong on the SSD, or at least on a drive with room. The project folder can be on an HDD:
  the work files are on Docker's disk, so only the results are written there (and `intermediate_files\`
  only with `publish_intermediates: true`).
- **Free space**: Docker data several GB (grows with each image pull); the image itself is
  several GB. Each dataset needs room on Docker's disk for its work files and room for the
  results on the project's drive (the final BAMs are the largest part; see README "Disk
  space"). The raw FASTQs stay where they are.
- **Virtualization**: `False` means it must be enabled in the BIOS/UEFI (Intel VT-x / AMD
  SVM). That is for the user or their IT, not you.
- **Memory hogs**: note big consumers (browsers with many tabs, IDE language servers,
  R or Python sessions with parallel workers). See section 5.

---

## 2. Windows

### 2.1 Install WSL 2 and Docker Desktop
Easiest: the pipeline does it. Clone the pipeline folder (section 3) and run
`RNASEQ\WINDOWS\run_rnaseq.bat setup`. If Docker Desktop is missing it
asks, then runs `wsl.exe --install --no-distribution` and
`winget install -e --id Docker.DockerDesktop`.

- **Administrator prompts (UAC)**: Windows shows them on the user's screen; tell the user
  to approve. A prompt that nobody approves makes the command hang or fail silently.
- **Restart required** after the first WSL/Docker install. Then the user opens **Docker
  Desktop** once from the Start menu, accepts the terms, and may skip sign-in ("Skip": no
  Docker account needed). Wait for "Engine running". Until that first launch, `docker info`
  fails and nothing works.
- **Where Docker keeps its data** is decided at install time. With winget it goes to `C:`
  (`%LOCALAPPDATA%\Docker\wsl`). If `C:` is small or not an SSD, install with an explicit
  data root instead (elevated, from the downloaded installer):
  `"Docker Desktop Installer.exe" install --accept-license --backend=wsl-2 --wsl-default-data-root=E:\DockerData`
  The setting ends up in `C:\ProgramData\DockerDesktop\install-settings.json`
  (`wslDefaultDataRoot`). Never pick a hard disk drive if an SSD has room. The work volume
  lives there too, so the drive needs room for the work files.
- If `winget` is missing (older Windows/LTSC): download Docker Desktop from docker.com and
  choose "Use WSL 2".

### 2.2 Check
```powershell
docker info --format '{{.NCPU}} CPUs, {{.MemTotal}} bytes for Docker'
docker run --rm hello-world        # needs internet; proves pulls and containers work
```
If `docker` is "not recognized" right after installing: the terminal predates the install.
Open a new one (the launcher adds Docker's bin folder itself).

### 2.3 Memory available to Docker (the launcher handles it; leave it alone)
Docker Desktop on WSL 2 takes its memory from `memory=` in `%USERPROFILE%\.wslconfig`
(`[wsl2]` section). By default it is half of the computer's memory.

For an RNA-seq run, the launcher sets this limit for the analysis, by the computer's total
memory: below 12 GB, 6 GB; 12 to under 24 GB, 8 GB; 24 GB or more, 16 GB. Then:

- it **asks first** (`--yes` accepts). Ask the user yourself before you run an analysis
  with `--yes`: the change restarts Docker Desktop and runs `wsl --shutdown`, which
  stops every WSL distribution on the computer;
- it changes **only** the `memory=` line and writes a marker comment
  (`# run_rnaseq: previous memory=<old or none>`); it keeps a copy of the file as
  `%USERPROFILE%\.wslconfig.run_rnaseq.bak` while the change is in place;
- it does this **only when no container is running** (`docker ps -q` empty); otherwise it
  keeps the current limit and says so;
- after the analysis (success or failure), if nothing else is running, it puts the old
  value back the same way and removes the marker and the copy. If the window was closed
  first, the next `run_rnaseq` command that finds Docker running and nothing else running
  restores it.
- `--no-wsl-resize` (on `run`) or `wsl_resize: false` in `config.yaml` leaves the file alone.

So do **not** raise Docker's memory "to be safe" by hand: Linux fills whatever the VM may use
with file cache, and Windows gets it back slowly. That, plus other programs, is what froze
the lab's first PC (see 2.5 and section 6). If the user has changed `.wslconfig` by hand,
the marker is absent and the launcher asks before it changes the value. Do not change a
value the user set by hand without their agreement.

### 2.4 Page file on an SSD
If `Win32_PageFileUsage` shows the page file on a hard disk drive, Windows pages to that slow
disk whenever memory runs short, and the whole computer stalls (worst when the same disk is
also busy with Docker or a FASTQ). With the user's OK, move it to the SSD as system-managed:
one elevated PowerShell run, effective after a restart:
```powershell
$cs = Get-CimInstance Win32_ComputerSystem
Set-CimInstance -InputObject $cs -Property @{ AutomaticManagedPagefile = $false }
Get-CimInstance Win32_PageFileSetting | Remove-CimInstance
New-CimInstance -ClassName Win32_PageFileSetting -Property @{ Name = 'C:\pagefile.sys'; InitialSize = [uint32]0; MaximumSize = [uint32]0 }
```
(Run it through `Start-Process powershell.exe -Verb RunAs -Wait -ArgumentList ...` so the
UAC prompt appears; check after the restart with `Get-CimInstance Win32_PageFileUsage`.)
Use the drive letter of the SSD in the `Name` if it is not `C:`.

### 2.5 Moving Docker's data later (lessons learned)
- Preferred: Docker Desktop > Settings > Resources > Advanced > **Disk image location** (the
  user clicks it; it moves the data itself).
- **Do not copy the `.vhdx` files by hand.** The copy loses a permission entry for WSL's
  virtual machine (a SID like `S-1-5-83-1-...`); Docker Desktop then fails with
  `Wsl/Service/AttachDisk/MountDisk/HCS/E_ACCESSDENIED` and an error dialog "disk not found:
  no sd* disk in /sys/block with wwid ...". Granting "NT VIRTUAL MACHINE\Virtual Machines"
  did not fix it.
- What worked: quit Docker Desktop (`Stop-Process` its processes) and run `wsl --shutdown`;
  set `wslDefaultDataRoot` in `C:\ProgramData\DockerDesktop\install-settings.json` to the new
  folder (elevated; keep a `.bak`); make sure no copied `disk\docker_data.vhdx` is in the new
  folder; start Docker Desktop: it creates a fresh, empty data disk there. Then run the
  pipeline's `setup` again to pull the image. This loses only images and caches: the lab
  pipelines keep nothing else inside Docker.
- The same applies to Docker's small system disk `main\ext4.vhdx`: WSL keeps its own
  registration of the `docker-desktop` distribution pointing at the OLD folder
  (`HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss\*`, `BasePath`). After changing the
  data root, quit Docker Desktop, run `wsl --shutdown` and `wsl --unregister docker-desktop`
  (it only holds Docker's own system files; the data disk with the images is separate), make
  sure the new folder contains no hand-copied `main\ext4.vhdx`, then start Docker Desktop: it
  re-registers the distribution in the new folder. Check with `wsl -l -v` and the `.vhdx`
  timestamps in the new folder, and **only then** delete the old folder.

### 2.6 Notes for the agent's own shells on Windows
- Agents typically have Git Bash and Windows PowerShell 5.1. The pipeline's `.bat` launcher
  runs PowerShell 5.1 with `-ExecutionPolicy Bypass`, so the user's execution policy doesn't
  matter.
- PowerShell 5.1 strips double quotes inside arguments passed to programs (`docker inspect
  --format '{{index .Labels "x"}}'` reaches docker without them). Use `{{json ...}}` plus
  `ConvertFrom-Json`, or Git Bash.
- `| Select-Object -First N` after a native command stops that command early (exit 255,
  half-done work). Don't pipe long-running launchers into it.
- Git Bash rewrites `/paths` in arguments to Windows paths; for `docker run` and `docker cp`
  use `MSYS_NO_PATHCONV=1` and convert host paths with `cygpath -w`. Git Bash's `grep` hides
  `\r`: count CRs with `tr -cd '\r' | wc -c`.
- `python` may be the Microsoft Store alias and hang. Use the pipeline image's Python
  (`docker run --rm ghcr.io/genesig/omics:<tag> python3 ...`) or Git Bash's if present.
- Excel files: if Excel is installed, convert with its COM object
  (`New-Object -ComObject Excel.Application`, `SaveAs(<csv>, 62)` = CSV UTF-8); otherwise
  convert them inside the pipeline image.
- Background jobs that poll a status file can be fooled by the previous run's file for a few
  seconds after `run --detach`; the launcher writes `running` immediately, but prefer the
  launcher's `status` command.

---

## 3. Clone the pipeline and pull its image

- Clone where the user keeps tools, e.g. `D:\PIPELINES\`. Not under OneDrive. The labmate
  gets only this folder with a sparse clone:
  `git clone --filter=blob:none --sparse <url>`, then `cd PIPELINES`, then
  `git sparse-checkout set RNASEQ/WINDOWS`. The repository root carries a
  `.gitattributes` that keeps Linux scripts LF even with Windows' `core.autocrlf=true`; it
  always comes along with a sparse checkout. On Windows also run `git config core.fileMode
  false` inside the clone (Windows has no executable bit).
- Run the pipeline's setup once: `RNASEQ\WINDOWS\run_rnaseq.bat setup --yes`.
  It **pulls** the pipeline's image, `ghcr.io/genesig/omics:<tag>`, with `docker pull
  --platform linux/amd64` (the tag is the one line of `docker\image.txt`). There is **no
  image build** on the user's computer: no Dockerfile in this folder. It takes several GB of
  download (keep the computer online and awake; corporate proxies and VPNs can break it), then
  prints the versions of the tools (Nextflow, HISAT2, samtools and the others) and adds the
  command to the user's PATH (new terminals only). Run it in the background or with a long
  timeout. `--yes` answers yes to routine questions. `setup` does not change `.wslconfig`.
- `run` pulls the image on its own if it is missing, so a failed pull in `setup` is fixed by
  running `setup` again.
- Validation: confirm that `docker image inspect ghcr.io/genesig/omics:<tag>` works and that
  the `setup` version lines print. Then check the pipeline on the benchmark (`README.md`,
  "Checking a run against the benchmark") before you trust a dataset's numbers.
- The image's tag changes only when `docker\image.txt` changes (`git pull`). A new tag is
  pulled by the next `setup` (or `run`).

---

## 4. Where data and projects go

Example layout on the lab's PC (the layout is a suggestion, not a requirement; the paths
are whatever the user writes in `config.yaml`):

```
D:\FASTQ\<run name>\                  raw sequencing runs: one folder per run, the FASTQ files
                                      (R1, R2) and run-level files; read in place, never moved
D:\REFERENCES\                        shared reference files, used by several datasets
  hisat2\t2t\                         HISAT2 index prefix of the genome (chm13-2.1.ht2 ... chm13-2.8.ht2)
  genomes\t2t\                        gene annotation GTF of the same build
D:\RNASEQ\<dataset>\                  analysis projects: one folder per dataset (made by `init`)
```

- Raw sequencing runs stay in the FASTQ folder. Several projects can point at the same
  run; the pipeline reads it in place and never moves, copies or renames it. Keep the path
  stable; if the run folder moves, update each project's `config.yaml` and samplesheet.
  Network shares (`\\server\...`) can't be read by Docker Desktop: copy the run to a local or
  USB drive first. The launcher refuses UNC paths and says so.
- Reference files are used in place too. `config.yaml` holds their paths (for example
  `hisat2_index: D:\REFERENCES\hisat2\t2t\chm13-2`).
- Analysis projects: one folder per dataset, created by the pipeline's `init` (e.g.
  `D:\RNASEQ\<dataset>\`), on an SSD if possible. If the chosen folder already exists and
  holds anything, list it and ask before using it (`init` adds files and a git repository to
  whatever is there). The project folder must be outside the pipeline folder; the launcher
  refuses otherwise.

---

## 5. Before each long run: free memory and disk

The pipeline caps its own memory, but it can't stop other programs. On the lab's first PC an
IDE's language server grew to 12.5 GB and R sessions held 32 worker processes (5 GB): with
Docker on top, Windows fell to 2 GB free and stalled. Before a long run (and if a run or
`docker` commands hang), check the top processes (section 1) and ask the user to close or
restart the big ones.

Before a run:
- check that no other container is running (`docker ps`); the launcher changes Docker's
  memory limit only when none is running, and otherwise runs with the limit Docker has;
- check the free disk space on Docker's disk (the drive with Docker Desktop's disk image,
  which holds the work volume) and on the project's drive against the launcher's figures
  (the launcher prints them at the start of the run; the project's drive stops below 20 GB free);
- check Docker's disk is on an SSD or a drive with room; the project folder may be on a hard
  disk (the launcher warns), since only results are written there with the work volume;
- confirm with the user that the Docker memory change and the Docker restart are acceptable
  (section 2.3), and that the computer will stay awake.

The launcher prints the budget at the start of each run (memory and cores). If the budget
is below the pipeline's minimum, it stops and says why: free memory, close programs, or use a
larger computer. Do not raise `max_memory_gb` above what the computer can spare.

---

## 6. Troubleshooting (symptom → cause → fix)

| Symptom | Likely cause | Fix |
|---|---|---|
| Whole PC freezes during an image pull or a run; event log: services time out, later "unexpected shutdown"; Docker VM log (`%LOCALAPPDATA%\Docker\log\vm\init.log`): "soft lockup", "page_reporting" | Memory pressure (Docker VM file cache + other programs) with the page file on a busy hard disk | Section 5 (close hogs), page file to SSD (2.4), Docker data on SSD (2.1/2.5) |
| `docker` commands hang for minutes; `status` hangs | Docker Desktop starved (low memory, disk saturated) | Free memory; wait; if needed quit Docker Desktop and `wsl --shutdown`, start it again |
| "Docker Desktop is unable to start"; `E_ACCESSDENIED` attaching `docker_data.vhdx` | Data disk copied by hand | 2.5 |
| "Docker Desktop is unable to start" after a memory change or after Docker sat idle; its log (`%LOCALAPPDATA%\Docker\log\host\com.docker.backend.exe.log`) shows `DockerDesktop/Wsl/CommandTimedOut` for `wsl.exe -l -v --all` | The WSL service hung while Docker was waking its VM | Quit Docker Desktop (end its processes), `wsl --shutdown`, start Docker Desktop again; then run the analysis again |
| Analysis says Docker's limit "stays at ... because other containers are running" | Another container (another pipeline, a test) is running | Stop it, then run again; the limit is changed only when nothing else runs |
| `%USERPROFILE%\.wslconfig.run_rnaseq.bak` exists and the `memory=` line is still the analysis value | The window was closed before the restore | Run any `run_rnaseq` command with Docker running and nothing else running; it puts the old value back |
| `run` refuses: budget below the pipeline's minimum | Too little memory for the computer, or `max_memory_gb` set low | Free memory, raise `max_memory_gb` if the computer can spare it, or use a larger computer |
| A step ends with exit 137 (killed for memory) | The step needed more than its budget (retried up to three times, never above the budget) | Free memory and run again; set `max_memory_gb` higher only if the computer can spare it |
| "names a network path" / UNC refused | Dataset on a network share | Copy to a local or USB drive, or map the share to a drive letter and use that path |
| "drive ... is not connected" | USB drive or mapped drive missing | Connect it, or correct the path in `config.yaml`/the samplesheet |
| `setup` cannot pull the image (network) | No internet, proxy, VPN | Fix connectivity; run `setup` again |
| Analysis much slower than expected; Windows sluggish | Docker's disk image on a hard disk drive, or the project on a hard disk with `work_volume: false` (the launcher warns for the project) | Move Docker's disk image to the SSD (Docker Desktop > Settings > Resources > Advanced > Disk image location); with `work_volume: false`, move the project too |

---

## 7. Hand-off to the user

Tell them, briefly: Docker Desktop must be running (it starts with the computer only if
they enable that in its settings); keep the computer plugged in and awake during analyses;
close big programs before long runs; the launcher may restart Docker Desktop for a run and
put its memory setting back afterwards (it asks first); raw data lives in the FASTQ folder
and is never moved by the pipeline; how to start a new dataset (`init`, see `README.md`); and
that an AI agent can do all of it for them.
