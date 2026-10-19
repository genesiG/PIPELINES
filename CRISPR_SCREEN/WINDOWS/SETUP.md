# SETUP.md — one-time setup of a labmate's Windows PC for the lab's Docker pipelines

Written for an AI agent (e.g. Claude Code) preparing a labmate's **Windows PC**
to run the lab's containerized pipelines (CRISPR_SCREEN today; RNA-seq,
CUT&RUN later). Do this **once per computer**. Afterwards each pipeline only
needs its own `run_crisprscreen setup` (the pipeline's own setup command, which
builds that pipeline's Docker image) and then runs from its runbook
(`AGENTS.md` in this folder). Nothing here is specific to one pipeline.

The user is usually a wet-lab scientist: explain each step in plain language,
ask before anything that needs an administrator prompt, a restart, or changes
system settings, and never delete their files.

Everything below comes from setting up the first lab PC (Windows 11, 2026-10):
the hiccups listed really happened.

---

## 0. End state (what "set up" means)

1. Docker Desktop installed, its terms accepted, the engine running
   (`docker info` works), on WSL 2.
2. Docker's own data (its disk image: images, build cache) on an **internal SSD**
   with ≥ 30 GB free.
3. The page file (virtual memory) on an SSD, not on a hard disk drive.
4. The pipeline repository cloned (not inside OneDrive-synced folders),
   its `setup` run once, and its synthetic test project passing.
5. A folder for sequencing data, e.g. `D:\FASTQ\<run>\`: one subfolder per
   sequencing run, raw FASTQs stay there; analysis projects point at them
   (they are never moved or copied).
6. The user knows: keep the computer plugged in and awake during runs, close
   memory-hungry programs, and how to start/check a run (the pipeline runbook).

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
- **RAM**: 8 GB or more is comfortable; 6 GB works (the CRISPR pipeline needs
  2.5 GB and runs step by step); less than 6 GB: tell the user it will be slow
  or may not fit.
- **Disks**: which drive is an SSD (`MediaType SSD`/NVMe) and which is a hard
  disk (`HDD`). Docker's data, the page file and (ideally) the analysis projects
  belong on the SSD. Projects on an HDD work but demultiplexing is several times
  slower and Windows may respond slowly meanwhile.
- **Free space**: Docker data ~5-20 GB (grows with each pipeline image and
  rebuild); each project needs roughly 1-2x the size of its raw FASTQ (per-sample
  FASTQs + `work/` while a run is going; `work/` is removed after a successful run
  unless `keep_work: true`), the raw FASTQ itself stays in the FASTQ folder.
- **Virtualization**: `False` means it must be enabled in the
  BIOS/UEFI (Intel VT-x / AMD SVM) — that is for the user or their IT, not you.
- **Memory hogs**: note big consumers (browsers with many tabs, IDE language
  servers, R/Python sessions with parallel workers). See section 5.

---

## 2. Windows

### 2.1 Install WSL 2 and Docker Desktop
Easiest: the pipeline does it — clone the repository (section 3) and run
`WINDOWS\run_crisprscreen.bat setup`. If Docker Desktop is missing it asks,
then runs `wsl --install --no-distribution` and
`winget install -e --id Docker.DockerDesktop`.

- **Administrator prompts (UAC)**: Windows shows them on the user's screen;
  tell the user to approve. A prompt that nobody approves makes the command
  hang or fail silently.
- **Restart required** after the first WSL/Docker install. Then the user opens
  **Docker Desktop** once from the Start menu, accepts the terms, and may skip
  sign-in ("Skip" — no Docker account needed). Wait for "Engine running".
  Until that first launch, `docker info` fails and nothing works.
- **Where Docker keeps its data** is decided at install time. With winget it
  goes to `C:` (`%LOCALAPPDATA%\Docker\wsl`). If `C:` is small or not an SSD,
  install with an explicit data root instead (elevated, from the downloaded
  installer):
  `"Docker Desktop Installer.exe" install --accept-license --backend=wsl-2 --wsl-default-data-root=E:\DockerData`
  The setting ends up in `C:\ProgramData\DockerDesktop\install-settings.json`
  (`wslDefaultDataRoot`). Never pick a hard disk drive if an SSD has room.
- If `winget` is missing (older Windows/LTSC): download Docker Desktop from
  docker.com, choose "Use WSL 2".

### 2.2 Check
```powershell
docker info --format '{{.NCPU}} CPUs, {{.MemTotal}} bytes for Docker'
docker run --rm hello-world        # needs internet; proves pulls and containers work
```
If `docker` is "not recognized" right after installing: the terminal predates
the install — open a new one (the pipeline launcher adds Docker's bin folder
itself).

### 2.3 Memory available to Docker (usually leave as is)
Docker Desktop on WSL 2 may use half of the RAM by default. The lab pipelines
cap themselves (CRISPR: 8 GB at most), so the default is enough on any machine
with ≥ 18 GB. On smaller machines the pipeline's `setup` explains and offers to
raise it (`run_crisprscreen setup --docker-memory`, which writes
`%USERPROFILE%\.wslconfig` and restarts WSL + Docker, stopping anything else
running in WSL). Ask the user first; `--yes` never does this.

Do **not** give Docker more memory "to be safe": Linux fills whatever the VM may
use with file cache, and Windows gets it back slowly. That, plus other programs,
is what froze the first lab PC (see 2.5 and section 6).

### 2.4 Page file on an SSD
If `Win32_PageFileUsage` shows the page file on a hard disk drive, Windows
pages to that slow disk whenever memory runs short, and the whole computer
stalls (worst when the same disk is also busy with Docker or a FASTQ). With the
user's OK, move it to the SSD as system-managed — one elevated PowerShell run,
effective after a restart:
```powershell
$cs = Get-CimInstance Win32_ComputerSystem
Set-CimInstance -InputObject $cs -Property @{ AutomaticManagedPagefile = $false }
Get-CimInstance Win32_PageFileSetting | Remove-CimInstance
New-CimInstance -ClassName Win32_PageFileSetting -Property @{ Name = 'C:\pagefile.sys'; InitialSize = [uint32]0; MaximumSize = [uint32]0 }
```
(Run it through `Start-Process powershell.exe -Verb RunAs -Wait -ArgumentList ...`
so the UAC prompt appears; check after the restart with
`Get-CimInstance Win32_PageFileUsage`.)

### 2.5 Moving Docker's data later (lessons learned)
- Preferred: Docker Desktop > Settings > Resources > Advanced > **Disk image
  location** (the user clicks it; it moves the data itself).
- **Do not copy the `.vhdx` files by hand.** The copy loses a permission entry
  for WSL's virtual machine (a SID like `S-1-5-83-1-...`); Docker Desktop then
  fails with `Wsl/Service/AttachDisk/MountDisk/HCS/E_ACCESSDENIED` and an error
  dialog "disk not found: no sd* disk in /sys/block with wwid ...". Granting
  "NT VIRTUAL MACHINE\Virtual Machines" did not fix it.
- What worked: quit Docker Desktop (`Stop-Process` its processes) and run
  `wsl --shutdown`; set `wslDefaultDataRoot` in
  `C:\ProgramData\DockerDesktop\install-settings.json` to the new folder
  (elevated; keep a `.bak`); make sure no copied `disk\docker_data.vhdx` is in
  the new folder; start Docker Desktop — it creates a fresh, empty data disk
  there. Then rebuild each pipeline's image (`setup`). This loses only images
  and caches: the lab pipelines keep nothing else inside Docker.
- The same applies to Docker's small system disk `main\ext4.vhdx`: WSL keeps
  its own registration of the `docker-desktop` distribution pointing at the
  OLD folder (`HKCU\Software\Microsoft\Windows\CurrentVersion\Lxss\*`,
  `BasePath`). After changing the data root, quit Docker Desktop, run
  `wsl --shutdown` and `wsl --unregister docker-desktop` (it only holds Docker's
  own system files, rebuilt from the installer's image; the data disk with the
  images is separate), make sure the new folder contains no hand-copied
  `main\ext4.vhdx`, then start Docker Desktop: it re-registers the distribution
  in the new folder. Check with `wsl -l -v` and the `.vhdx` timestamps in the
  new folder, and **only then** delete the old folder. (On the first lab PC the
  old folder was deleted while the registration still pointed at it: Docker
  failed with `Wsl/Service/CreateInstance/MountDisk/ERROR_PATH_NOT_FOUND`, then,
  with a copied `ext4.vhdx`, `RegisterDistro/MountDisk/HCS/E_ACCESSDENIED`.)

### 2.6 Notes for the agent's own shells on Windows
- Agents typically have Git Bash and Windows PowerShell 5.1. The pipeline's
  `.bat` launcher runs PowerShell 5.1 with `-ExecutionPolicy Bypass`, so the
  user's execution policy doesn't matter.
- PowerShell 5.1 strips double quotes inside arguments passed to programs
  (`docker inspect --format '{{index .Labels "x"}}'` reaches docker without
  them). Use `{{json ...}}` + `ConvertFrom-Json`, or Git Bash.
- `| Select-Object -First N` after a native command stops that command early
  (exit 255, half-done work). Don't pipe long-running launchers into it.
- Git Bash rewrites `/paths` in arguments to Windows paths; for `docker run`
  and `docker cp` use `MSYS_NO_PATHCONV=1` and convert host paths with
  `cygpath -w`. Git Bash's `grep` hides `\r`: count CRs with `tr -cd '\r' | wc -c`.
- `python` may be the Microsoft Store alias and hang. Use the pipeline image's
  Python (`docker run --rm <image> python3 ...`) or Git Bash's if present.
- Excel files: if Excel is installed, convert with its COM object
  (`New-Object -ComObject Excel.Application`, `SaveAs(<csv>, 62)` = CSV UTF-8);
  otherwise R's readxl inside the pipeline image.
- Background jobs that poll a status file can be fooled by the previous run's
  file for a few seconds after `run --detach`; the launchers now write
  `running` immediately, but prefer the launcher's `status` command.

---

## 3. Clone the pipeline(s) and build each image

- Clone where the user keeps tools, e.g. `D:\PIPELINES\`. Not under OneDrive.
  The labmate obtained this folder with a sparse clone (only `WINDOWS/` is
  checked out):
  `git clone --filter=blob:none --sparse <url>`, then `cd CRISPR_SCREEN`, then
  `git sparse-checkout set WINDOWS`. The repository root carries a
  `.gitattributes` that keeps Linux scripts LF even with Windows'
  `core.autocrlf=true`; it always comes along with a sparse checkout. On Windows
  also run `git config core.fileMode false` inside the clone (Windows has no
  executable bit).
- Run the pipeline's `setup --yes` once (CRISPR: `WINDOWS\run_crisprscreen.bat
  setup --yes`): it builds the
  image (10-30 min, ~2-4 GB download, needs internet; corporate proxies/VPNs
  can break it), prints tool versions and adds the command to the user's PATH
  (new terminals only). Run it in the background or with a long timeout.
  `--yes` answers yes to routine questions but never edits `.wslconfig`.
- A transient micromamba crash during the build
  (`Unexpected error 9 on netlink descriptor`) is retried automatically; if the
  whole build fails, simply run `setup` again (finished layers are cached).
- Validate with the pipeline's synthetic test (`test\README.md` in this
  folder): a few minutes; it must end "succeeded" with the
  expected hits.
- Keep only the pipeline images: helper images pulled for maintenance
  (e.g. BioContainers used to regenerate lock files) can be removed with
  `docker rmi`.

## 4. Where data and projects go

- Raw sequencing runs: `D:\FASTQ\<run name>\` or an internal-SSD folder of the
  user's choice — one subfolder per run, with
  the multiplexed `Undetermined*_R1_001.fastq.gz` (and R2), plus run-level
  files (full sample sheet, run QC). Several labmates' projects can point to the
  same run; the pipelines read it in place and never move, copy or rename it.
  Keep the path stable; if the run folder moves, update each project's
  `config.yaml`. Network shares (`\\server\...`) can't be read by Docker
  Desktop: copy the run to a local or USB drive first.
- Analysis projects: one folder per screen, created by the pipeline's `init`
  (e.g. `D:\<SCREEN_NAME>\`), on an SSD if
  possible. If the chosen folder already exists and holds anything, list it and
  ask before using it (`init` adds files and a git repository to whatever is
  there).

## 5. Before each long run: free memory

The pipelines cap their own memory, but they can't stop other programs. On the
first lab PC an IDE's language server grew to 12.5 GB and R sessions held 32
worker processes (5 GB): with Docker on top, Windows fell to 2 GB free and
stalled. Before a long run (and if a run or `docker` commands hang), check the
top processes (section 1) and ask the user to close or restart the big ones.
The CRISPR launcher warns when free memory is below the run's budget, and on
Windows it makes Docker drop its file cache when Windows runs low.

## 6. Troubleshooting (symptom → cause → fix)

| Symptom | Likely cause | Fix |
|---|---|---|
| Whole PC freezes during an image build or a run; event log: services time out, later "unexpected shutdown"; Docker VM log (`%LOCALAPPDATA%\Docker\log\vm\init.log`): "soft lockup", "page_reporting" | Memory pressure (Docker VM file cache + other programs) with the page file on a busy hard disk | Section 5 (close hogs), page file to SSD (2.4), Docker data on SSD (2.1/2.5); the pipeline's memory cap and cache guard |
| `docker` commands hang for minutes; `status` hangs | Docker Desktop starved (low memory, disk saturated) | Free memory; wait; if needed quit Docker Desktop and `wsl --shutdown`, start it again |
| "Docker Desktop is unable to start"; `E_ACCESSDENIED` attaching `docker_data.vhdx` | Data disk copied by hand | 2.5 |
| "Docker Desktop is unable to start" after Docker sat idle for a while; its log (`%LOCALAPPDATA%\Docker\log\host\com.docker.backend.exe.log`) shows `DockerDesktop/Wsl/CommandTimedOut` for `wsl.exe -l -v --all` | The WSL service hung while Docker was waking its idle VM | Quit Docker Desktop (end its processes), `wsl --shutdown`, start Docker Desktop again (seconds); a run started meanwhile fails after the launcher's 4-minute wait — just start it again |
| `docker info` fails after install | First launch/terms not done, or restart pending | Restart, open Docker Desktop, accept terms |
| "Virtualization is disabled" / WSL errors | VT-x/SVM off in firmware | User/IT enables it in BIOS/UEFI |
| Image build fails (network) | No internet, proxy, VPN | Fix connectivity; rerun `setup` |
| The launcher rebuilds the image on every run | (fixed) PowerShell 5.1 quote stripping broke the image check | Update the pipeline repository |
| Analysis much slower than expected; Windows sluggish | Project or Docker data on a hard disk drive (the launcher warns) | Move the project/data to the SSD |

## 7. Hand-off to the user

Tell them, briefly: Docker Desktop must be running (it starts with the
computer only if they enable that in its settings); keep the computer plugged
in and awake during analyses; close big programs before long runs; raw data
lives in the FASTQ folder and is never moved by the pipelines; how to start a
new screen (the pipeline runbook's `init`), and that an AI agent can do all of
it for them.
