# run_atacseq.ps1 - ATAC-seq pipeline launcher, Docker version for Windows.
#
# Labmates (and AI agents) call it through run_atacseq.bat, which only
# exists so PowerShell's execution policy never gets in the way:
#
#   run_atacseq setup            one-time: check/install Docker Desktop, download the
#                                  pipeline's Docker image, put this command on PATH
#   run_atacseq init [FOLDER]    set up a new project folder (default: the current folder)
#   run_atacseq [run] [CONFIG]   run the pipeline for the project in the current folder
#                                  (CONFIG defaults to config.yaml)
#   run_atacseq status           is it running? how far along? did it fail, and where?
#   run_atacseq logs             follow the live log of the running analysis
#   run_atacseq stop             stop a running analysis (re-running resumes it)
#   run_atacseq clean            delete the work folder, resume cache and Docker work volume now
#                                  (metadata\work, metadata\.nextflow, and the volume; results\ is kept)
#   run_atacseq script utils\X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment
#   run_atacseq help
#
# Options: --detach (run: start and return, don't follow the log), --yes
# (answer yes to questions: for AI agents / scripts), --pause (wait for Enter
# before closing; used by the double-clickable RUN_ANALYSIS.bat), --no-wsl-resize
# (run: leave Docker Desktop's memory limit alone; see Invoke-WslStep).
#
# What `run` does (the Docker counterpart of the cluster's run_atacseq.sh):
#   1. Checks the input files named in config.yaml and samplesheet.csv. The inputs
#      stay where they are: never moved, copied or renamed. A Windows drive they
#      name (D:\FASTQ\x.fastq.gz) is mounted read-only in the container. Network
#      paths (\\server\share) are refused: map the share to a drive letter first.
#   2. Checks the free disk space on the project's drive (warns; stops only below
#      disk_hard_min_gb in docker\requirements.conf).
#   3. Sets Docker Desktop's memory limit for this analysis, sized for this
#      computer (asks first; the old value is put back when the analysis is done).
#   4. Makes sure the pipeline's image is downloaded, then starts one container
#      that runs Nextflow and every step (docker/run_in_container.sh), with the
#      project folder mounted at /project and this pipeline folder (read-only) at
#      /pipeline, then follows its log. Closing the window does NOT stop the analysis.
#      The work folder (every step's working files and the resume cache) is a Docker
#      volume by default (work_volume: true in config.yaml; false keeps metadata\work
#      in the project folder). Only results\ is written to the project folder then.
#
# Hidden commands (maintainers and tests): wsl-plan, scan-paths, selftest.
# -WhatIfWsl (with run) prints the Docker memory plan and stops; it changes nothing.

$ErrorActionPreference = 'Continue'
Set-StrictMode -Version 2

$PipelineDir  = $PSScriptRoot
$PipelineName = 'atacseq'
$ImageFile    = Join-Path $PipelineDir 'docker\image.txt'
$Image        = ''                      # read from docker\image.txt by the commands that need it
$script:WslFile   = Join-Path $env:USERPROFILE '.wslconfig'
$script:WslBackup = Join-Path $env:USERPROFILE '.wslconfig.run_atacseq.bak'
$script:WslRestorePending = $false

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
$script:Pause = $false; $script:Yes = $false; $script:NoWslResize = $false; $script:WhatIfWsl = $false
$Detach = $false; $SelfTest = $false
$positional = @()
$rawRest = $false
foreach ($a in $args) {
    # After 'script' every argument goes to the script verbatim (even --yes or --pause).
    if ($rawRest) { $positional += [string]$a; continue }
    switch -Regex ([string]$a) {
        '^--?pause$'           { $script:Pause = $true; continue }
        '^--?yes$|^-y$'        { $script:Yes = $true; continue }
        '^--?detach$'          { $Detach = $true; continue }
        '^--?no-?wsl-?resize$' { $script:NoWslResize = $true; continue }
        '^-WhatIfWsl$'         { $script:WhatIfWsl = $true; continue }
        '^--?selftest$'        { $SelfTest = $true; continue }
        default                { $positional += [string]$a }
    }
    if ($positional.Count -eq 1 -and $positional[0] -eq 'script') { $rawRest = $true }
}
$Command = 'run'
if ($SelfTest) {
    $Command = 'selftest'
} elseif ($positional.Count -gt 0) {
    if ($positional[0] -match '^(setup|init|run|status|logs|stop|clean|script|help|wsl-plan|scan-paths|selftest|-h|--help|/\?)$') {
        $Command = $positional[0]
        $positional = @($positional | Select-Object -Skip 1)
    }
}

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
function Say([string]$msg)  { Write-Host $msg }
function Info([string]$msg) { Write-Host $msg -ForegroundColor Cyan }
function Warn([string]$msg) { Write-Host "WARNING: $msg" -ForegroundColor Yellow }
function Done([int]$code) {
    # Docker's memory limit is put back before the window closes (when nothing else is running).
    if ($script:WslRestorePending) { $script:WslRestorePending = $false; Restore-WslIfIdle }
    if ($script:Pause) { Write-Host ''; Read-Host 'Press Enter to close this window' | Out-Null }
    exit $code
}
function Fail([string]$msg) {
    Write-Host ''
    Write-Host "ERROR: $msg" -ForegroundColor Red
    Done 1
}
function Confirm([string]$question) {
    if ($script:Yes) { return $true }
    try { $ans = Read-Host "$question [y/N]" } catch { return $false }
    return ($ans -match '^(y|yes)$')
}
function Have([string]$cmd) { return [bool](Get-Command $cmd -ErrorAction SilentlyContinue) }

# Run a native command quietly; returns its exit code.
function Quiet([scriptblock]$block) {
    $null = & $block 2>&1
    return $LASTEXITCODE
}

function Write-TextLF([string]$path, [string]$text) {
    $utf8 = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($path, $text, $utf8)
}

# Byte-exact text I/O (no BOM added or removed; used for .wslconfig).
function Read-TextRaw([string]$path) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { return $null }
    $bytes = [System.IO.File]::ReadAllBytes($path)
    return (New-Object System.Text.UTF8Encoding($false)).GetString($bytes)
}
function Write-TextRaw([string]$path, [string]$text) {
    $bytes = (New-Object System.Text.UTF8Encoding($false)).GetBytes($text)
    [System.IO.File]::WriteAllBytes($path, $bytes)
}

function Format-GB([int]$mb) { return ('{0:0.#} GB' -f ($mb / 1024.0)) }
function Format-Bytes([double]$b) {
    if ($b -ge 1GB) { return ('{0:0.#} GB' -f ($b / 1GB)) }
    return ('{0:0} MB' -f ($b / 1MB))
}

# ---------------------------------------------------------------------------
# Docker
# ---------------------------------------------------------------------------
function Find-DockerDesktopExe {
    $candidates = @("$env:ProgramFiles\Docker\Docker\Docker Desktop.exe")
    try {
        $reg = Get-ItemProperty 'HKLM:\SOFTWARE\Docker Inc.\Docker\1.0' -ErrorAction Stop
        if ($reg.AppPath) { $candidates = @((Join-Path $reg.AppPath 'Docker Desktop.exe')) + $candidates }
    } catch { }
    foreach ($c in $candidates) { if (Test-Path $c) { return $c } }
    return $null
}

function Install-DockerDesktop {
    Say ''
    Say 'Docker Desktop is not installed. The pipeline runs inside Docker, so it is needed (one time).'
    if (-not (Have 'winget')) {
        Fail ("Please install Docker Desktop from https://www.docker.com/products/docker-desktop/ " +
              "(keep the default settings), restart the computer if asked, start Docker Desktop once, then run this again.")
    }
    if (-not (Confirm 'Install Docker Desktop now with winget? (Windows will ask for administrator permission)')) {
        Fail 'Docker Desktop is required. Install it from https://www.docker.com/products/docker-desktop/ and run this again.'
    }
    Info 'Installing Docker Desktop and the Windows parts it needs (this takes a while)...'
    & wsl.exe --install --no-distribution
    & winget install -e --id Docker.DockerDesktop --accept-package-agreements --accept-source-agreements
    if ($LASTEXITCODE -ne 0) { Fail "winget could not install Docker Desktop (exit $LASTEXITCODE)." }
    Say ''
    Say 'Docker Desktop is installed. Next:'
    Say '  1. RESTART the computer (needed once).'
    Say '  2. Open Docker Desktop from the Start menu and accept its terms (no account needed: "Skip").'
    Say '  3. Run this command again.'
    Done 0
}

# Makes sure the docker command exists and the Docker engine is running.
function Ensure-Docker([switch]$AllowInstall) {
    if (-not (Have 'docker')) {
        $exe = Find-DockerDesktopExe
        if ($exe) {
            # Installed, but this terminal's PATH predates it.
            $env:Path = "$env:ProgramFiles\Docker\Docker\resources\bin;$env:Path"
        }
    }
    if (-not (Have 'docker')) {
        if ($AllowInstall) { Install-DockerDesktop }
        Fail "Docker is not installed. Run:  run_atacseq setup"
    }
    if ((Quiet { docker info }) -eq 0) { return }

    $exe = Find-DockerDesktopExe
    if (-not $exe) { Fail 'Docker is installed but its engine is not running. Start Docker Desktop and try again.' }
    Info 'Starting Docker Desktop (this can take a minute)...'
    Start-Process -FilePath $exe | Out-Null
    $deadline = (Get-Date).AddMinutes(4)
    while ((Get-Date) -lt $deadline) {
        Start-Sleep -Seconds 5
        if ((Quiet { docker info }) -eq 0) { Say 'Docker is running.'; return }
    }
    Fail ('Docker Desktop did not start within 4 minutes. Open Docker Desktop yourself, wait until it says ' +
          '"Engine running" (accept its terms the first time), then run this again.')
}

function Test-DockerUp {
    if (-not (Have 'docker')) { return $false }
    return ((Quiet { docker info }) -eq 0)
}

# Containers running right now (any project).
function Get-RunningCount {
    if (-not (Have 'docker')) { return 0 }
    $ids = @(& docker ps -q 2>$null | Where-Object { $_ })
    return $ids.Count
}

# The image tag: the one line of docker\image.txt (ATAC_IMAGE_OVERRIDE replaces it).
function Get-ImageTag {
    if ($env:ATAC_IMAGE_OVERRIDE) { return $env:ATAC_IMAGE_OVERRIDE.Trim() }
    if (-not (Test-Path -LiteralPath $ImageFile)) {
        Fail "docker\image.txt is missing from the pipeline folder ($PipelineDir). It names the pipeline's software image."
    }
    foreach ($line in [System.IO.File]::ReadAllLines($ImageFile)) {
        $t = $line.Trim()
        if ($t -and -not $t.StartsWith('#')) { return $t }
    }
    Fail 'docker\image.txt is empty. It must hold the image name, e.g. ghcr.io/genesig/omics:2026.10.2'
}

# Downloads the pipeline's image (setup always does; run only when it is missing).
function Pull-Image {
    Info "Downloading the pipeline's software image ($Image)."
    Say  'This happens once (and again only when docker\image.txt names a new version).'
    Say  'It is several GB and can take a while. Keep the computer awake and online.'
    Start-Guard -watchPid $PID
    & docker pull --platform linux/amd64 $Image
    if ($LASTEXITCODE -ne 0) {
        Fail ("Downloading the pipeline's image failed (see the messages above). Check the internet connection " +
              "and free disk space, then run this again.")
    }
    Info "Image $Image is ready."
}

function Ensure-Image {
    if ((Quiet { docker image inspect $Image }) -eq 0) { return }
    Pull-Image
}

# ---------------------------------------------------------------------------
# Pipeline files: bash scripts must have LF line endings inside the Linux
# container. .gitattributes keeps them LF on checkout; this repairs a copy
# that was made some other way (zip download, copy through a CRLF tool).
# ---------------------------------------------------------------------------
function Repair-LineEndings {
    $targets = @()
    $targets += @(Get-ChildItem (Join-Path $PipelineDir 'bin') -Filter '*.sh' -ErrorAction SilentlyContinue)
    $targets += @(Get-ChildItem (Join-Path $PipelineDir 'docker') -Filter '*.sh' -Recurse -ErrorAction SilentlyContinue)
    foreach ($f in $targets) {
        $text = [System.IO.File]::ReadAllText($f.FullName)
        if ($text.Contains("`r")) {
            Write-TextLF $f.FullName ($text -replace "`r`n", "`n" -replace "`r", "`n")
            Warn "converted $($f.Name) to Unix line endings (needed inside the container)"
        }
    }
}

# ---------------------------------------------------------------------------
# config.yaml: flat keys only (as on the cluster), so simple line reads do.
# ---------------------------------------------------------------------------
function Get-YamlValue([string]$file, [string]$key) {
    if (-not (Test-Path -LiteralPath $file)) { return $null }
    foreach ($line in [System.IO.File]::ReadAllLines($file)) {
        if ($line -match "^$([regex]::Escape($key)):(.*)$") {
            $v = $Matches[1] -replace '\s+#.*$', ''
            $v = $v.Trim()
            if ($v -match '^"(.*)"$' -or $v -match "^'(.*)'$") { $v = $Matches[1] }
            return $v
        }
    }
    return $null
}

function Test-Unset([string]$v) { return (-not $v) -or ($v -eq 'null') -or ($v -eq '~') }

function Get-ConfigBool([string]$key, [bool]$default) {
    $v = Get-YamlValue $script:ConfigFile $key
    if (Test-Unset $v) { return $default }
    if ($v -match '^(false|no|off|0)$') { return $false }
    return $true
}

# ---------------------------------------------------------------------------
# Docker work volume (work_volume: true, the default): Nextflow's work folder and
# resume cache live in a named Docker volume (Docker's own disk, much faster than
# the project folder's bind mount), not in metadata\. Pure helpers first.
# ---------------------------------------------------------------------------
# One volume per project. Lower-cased first, so the same folder in another case
# always gets the same name.
function Get-WorkVolumeName([string]$dir) {
    return (Get-ContainerName ("$dir".ToLowerInvariant())) + '_work'
}

# '' when the publish mode is fine with the work volume, else the message to show.
function Get-PublishModeProblem([string]$mode, [bool]$workVolume) {
    if (-not $workVolume) { return '' }
    $m = ("$mode").Trim().ToLowerInvariant()
    if (@('symlink', 'rellink', 'link') -contains $m) {
        return ("publish_mode: $m cannot be used with the Docker work volume. Hard links and links cannot " +
                "cross filesystems (the work folder is on Docker's disk, results\ is in the project folder), " +
                "and results\ would point into the volume. Set publish_mode: auto or copy in config.yaml, " +
                "or set work_volume: false.")
    }
    return ''
}

# Free bytes from `df -Pk <path>` text: the Available column (4th) of the last line.
function Get-DfAvailBytes([string]$dfText) {
    $lines = @($dfText -split "`r?`n" | Where-Object { $_.Trim() -ne '' })
    if ($lines.Count -eq 0) { return $null }
    $parts = @($lines[-1].Trim() -split '\s+')
    if ($parts.Count -lt 4 -or $parts[3] -notmatch '^\d+$') { return $null }
    return ([double]$parts[3]) * 1024
}

# How the work files fit in the free space: 'fail' when they do not fit, 'warn' when they fit with
# less than 20% to spare, '' when there is room. Pure (no disk access).
function Get-DiskFit([double]$freeB, [double]$needB) {
    if ($freeB -lt $needB) { return 'fail' }
    if ($freeB -lt ($needB * 1.2)) { return 'warn' }
    return ''
}

# The disk-image folder ("DataFolder") from Docker Desktop's settings JSON, or $null.
function Get-DockerDataFolder([string]$settingsJson) {
    if (-not $settingsJson) { return $null }
    try { $obj = $settingsJson | ConvertFrom-Json } catch { return $null }
    if ($null -eq $obj) { return $null }
    foreach ($p in @($obj.PSObject.Properties)) {
        if ($p.Name -ieq 'DataFolder' -and $p.Value -is [string] -and $p.Value.Trim() -ne '') { return $p.Value.Trim() }
    }
    return $null
}

# ---------------------------------------------------------------------------
# Input paths. Inputs are never moved: a Windows drive named in config.yaml or
# in the project's CSV files is mounted read-only at /mnt/<letter>, and
# main.nf reads D:\x\y as /mnt/d/x/y. Network paths cannot be mounted.
# ---------------------------------------------------------------------------
# A drive path starts at the start of a cell/value (or after a quote, comma,
# space or tab). Group 1 is the drive letter.
$script:DrivePattern = '(?:^|(?<=[\s,''"]))([A-Za-z]):[\\/][^\s,''"]*'
$script:UncPattern   = '(?:^|(?<=[\s,''"]))((?:\\\\|//)[^\\/\s,''"]+[\\/][^\s,''"]*)'

# A path as a full Windows path (relative = relative to the project folder).
function Resolve-InputPath([string]$v) {
    if ($v -match '^[A-Za-z]:[\\/]') { return [System.IO.Path]::GetFullPath(($v -replace '/', '\')) }
    return [System.IO.Path]::GetFullPath((Join-Path $script:ProjectDir ($v -replace '/', '\')))
}

# Every drive path and network path in config.yaml (comments removed) and in the
# project's CSV files (root folder only). Returns hashtables: File, Path, Kind, Drive.
function Get-PathRefs([string]$file) {
    $refs = New-Object System.Collections.ArrayList
    $isCsv = ($file -like '*.csv')
    foreach ($line in [System.IO.File]::ReadAllLines($file)) {
        $text = $line
        if (-not $isCsv) { $text = $line -replace '(^|\s)#.*$', '' }
        foreach ($m in [regex]::Matches($text, $script:UncPattern)) {
            [void]$refs.Add(@{ File = $file; Path = $m.Groups[1].Value; Kind = 'unc'; Drive = '' })
        }
        foreach ($m in [regex]::Matches($text, $script:DrivePattern)) {
            [void]$refs.Add(@{ File = $file; Path = $m.Value; Kind = 'drive'; Drive = $m.Groups[1].Value.ToUpperInvariant() })
        }
    }
    return $refs.ToArray()
}

function Get-ProjectPathRefs {
    $files = @()
    if (Test-Path -LiteralPath $script:ConfigFile) { $files += $script:ConfigFile }
    $files += @(Get-ChildItem -LiteralPath $script:ProjectDir -Filter '*.csv' -File -ErrorAction SilentlyContinue | ForEach-Object { $_.FullName })
    $all = @()
    foreach ($f in $files) { $all += @(Get-PathRefs $f) }
    return $all
}

# Stops on network paths and on drives that are not connected. Returns the drive letters to mount.
function Check-InputDrives {
    $refs = @(Get-ProjectPathRefs)
    $unc = @($refs | Where-Object { $_.Kind -eq 'unc' })
    if ($unc.Count -gt 0) {
        $u = $unc[0]
        Fail ("$(Split-Path $u.File -Leaf) names a network path ($($u.Path)). Docker Desktop cannot read network shares. " +
              "Map the share to a drive letter (File Explorer: This PC > Map network drive), then write the path with that letter, " +
              "for example Z:\data\x.fastq.gz.")
    }
    $drives = @{}
    foreach ($r in @($refs | Where-Object { $_.Kind -eq 'drive' })) {
        if (-not $drives.ContainsKey($r.Drive)) { $drives[$r.Drive] = $r }
    }
    $letters = @()
    foreach ($letter in @($drives.Keys | Sort-Object)) {
        $r = $drives[$letter]
        if (-not (Test-Path -LiteralPath "$($letter):\")) {
            Fail ("$(Split-Path $r.File -Leaf) names $($r.Path), but the drive $($letter): is not connected. " +
                  "Connect the drive (or correct the path) and run this again.")
        }
        $letters += $letter
    }
    return $letters
}

# docker run arguments that mount each drive read-only, and the ATAC_MOUNTS_HOST text (lines like D:=/mnt/d).
function Get-MountArgs([string[]]$letters) {
    $a = @()
    foreach ($l in $letters) { $lc = $l.ToLowerInvariant(); $a += @('-v', "${l}:\:/mnt/${lc}:ro") }
    return $a
}
function Get-MountsHostText([string[]]$letters) {
    $lines = @()
    foreach ($l in $letters) { $lc = $l.ToLowerInvariant(); $lines += "${l}:=/mnt/${lc}" }
    return ($lines -join "`n")
}

# Plain warnings only: the pipeline itself checks every file.
function Check-ConfigFiles {
    $keys = @('genome_index', 'contamination_index', 'blacklist', 'chrom_sizes', 'gtf', 'tss_bed')
    foreach ($k in $keys) {
        $v = Get-YamlValue $script:ConfigFile $k
        if (Test-Unset $v) { continue }
        $items = @($v)
        foreach ($item in $items) {
            if ($item -match '^[/~]') {
                Warn ("${k}: $item is not a Windows path. Write a path like D:\ref\x.bed, or a path relative to the project folder. " +
                      "The pipeline would look for it inside Docker and probably fail.")
                continue
            }
            $full = Resolve-InputPath $item
            if ($k -eq 'genome_index' -or $k -eq 'contamination_index') {
                $dir = Split-Path $full -Parent
                $leaf = Split-Path $full -Leaf
                $n = 0
                if (Test-Path -LiteralPath $dir -PathType Container) {
                    $n = @(Get-ChildItem -LiteralPath $dir -Filter "$leaf*.bt2*" -File -ErrorAction SilentlyContinue).Count
                }
                if ($n -eq 0) { Warn "${k}: no $leaf*.bt2* index files found in $dir" }
            } elseif (-not (Test-Path -LiteralPath $full -PathType Leaf)) {
                Warn "${k}: file not found: $item (looked for $full)"
            }
        }
    }
}

# The output folders are written inside the project folder; a Windows drive path would be read-only.
function Check-OutputDirs {
    $defaults = @{ outdir = './results'; metadata_dir = './metadata'; intermediate_dir = './intermediate_files' }
    foreach ($k in @('outdir', 'metadata_dir', 'intermediate_dir')) {
        $v = Get-YamlValue $script:ConfigFile $k
        if ((-not (Test-Unset $v)) -and ($v -match '^[A-Za-z]:[\\/]')) {
            Fail ("${k}: $v is a Windows path. The results are written into the project folder, so write it relative " +
                  "to the project folder instead, e.g. ${k}: $($defaults[$k]).")
        }
    }
}

function Get-SamplesheetFile {
    $v = Get-YamlValue $script:ConfigFile 'samplesheet'
    if (Test-Unset $v) { $v = 'samplesheet.csv' }
    return (Resolve-InputPath $v)
}

function Check-Samplesheet {
    $ss = Get-SamplesheetFile
    if (-not (Test-Path -LiteralPath $ss -PathType Leaf)) {
        Fail ("the samplesheet ($ss) is not in the project folder. Copy the template from examples\samplesheet.csv " +
              "into the project folder, fill it in, and run this again.")
    }
}

# Free space on the project's drive against the input files named in the samplesheet.
# The drive folder that holds Docker Desktop's disk image, as its settings name it (or $null).
function Get-DockerDiskFolder {
    foreach ($f in @('settings-store.json', 'settings.json')) {
        $p = Join-Path $env:APPDATA "Docker\$f"
        if (-not (Test-Path -LiteralPath $p -PathType Leaf)) { continue }
        $d = Get-DockerDataFolder ([System.IO.File]::ReadAllText($p))
        if ($d) { return $d }
    }
    return $null
}

# The drive that holds Docker's disk image: the drive of Docker's DataFolder when its settings name
# one, else the drive of %LOCALAPPDATA% (Docker Desktop's default). '' when neither is on a drive letter.
function Get-DockerDiskRoot {
    try {
        $dataFolder = Get-DockerDiskFolder
        if ($dataFolder) {
            $r = [System.IO.Path]::GetPathRoot($dataFolder)
            if ($r -match '^[A-Za-z]:\\$') { return $r }
        }
        if ($env:LOCALAPPDATA) {
            $r = [System.IO.Path]::GetPathRoot($env:LOCALAPPDATA)
            if ($r -match '^[A-Za-z]:\\$') { return $r }
        }
    } catch { }
    return ''
}

# Free space on Docker's disk, where the work volume lives. It is the free space of the drive that holds
# Docker's disk image (Get-DockerDiskRoot), not the df inside the image: that reports the VM's virtual disk,
# which can be far bigger than the drive. The df value is used only when the drive is unknown.
# The probing is silent on any failure. The decision (Fail or Warn) is made outside the try block, so a
# Fail is not swallowed.
function Check-DockerDisk([double]$needB, [bool]$onVolume) {
    if ((-not $onVolume) -or $script:WhatIfWsl) { return }
    $free = $null
    $root = ''
    try {
        if (-not (Test-DockerUp)) { return }
        if ((Quiet { docker image inspect $Image }) -ne 0) { return }
        $dfOut = @(& docker run --rm --platform linux/amd64 --memory 256m $Image df -Pk /tmp 2>$null)
        $freeDocker = Get-DfAvailBytes ($dfOut -join "`n")
        $root = Get-DockerDiskRoot
        $driveFree = $null
        if ($root) { $driveFree = [double](New-Object System.IO.DriveInfo -ArgumentList $root).AvailableFreeSpace }
        if ($null -ne $driveFree) { $free = $driveFree } else { $free = $freeDocker }
    } catch { return }
    if ($null -eq $free) { return }
    $drive = ''
    if ($root) { $drive = $root.TrimEnd('\') }
    $where = ''
    if ($drive) { $where = " (on drive $drive)" }
    Say "  Docker's disk: $(Format-Bytes $free) free$where"
    $fit = Get-DiskFit $free $needB
    if ($fit -eq 'fail') {
        $hostText = "the drive that holds Docker's disk"
        if ($drive) { $hostText = "drive $drive" }
        Fail ("Docker's disk has $(Format-Bytes $free) free$where, but the work files need about $(Format-Bytes $needB). " +
              "Set work_volume: false in $($script:ConfigName) (the work folder then goes into the project folder, on a drive " +
              "with more room), or free space on $hostText (run_atacseq clean deletes the Docker work volume), then run this again.")
    }
    if ($fit -eq 'warn') {
        Warn ("Docker's disk has $(Format-Bytes $free) free$where, and the work files may need about $(Format-Bytes $needB), " +
              "so there is little room to spare. If it fills up the run stops and can be resumed after you make room. " +
              "To free space: run_atacseq clean deletes the Docker work volume; keep_work: false in config.yaml deletes it " +
              "by itself after a successful run.")
    }
}

function Check-Disk {
    $ss = Get-SamplesheetFile
    $seen = @{}
    $total = [double]0
    $count = 0
    foreach ($row in @(Import-Csv -LiteralPath $ss)) {
        foreach ($col in @('fastq_1', 'fastq_2')) {
            if (-not ($row.PSObject.Properties.Name -contains $col)) { continue }
            $v = "$($row.$col)".Trim()
            if (Test-Unset $v) { continue }
            if ($v -match '^(\\\\|//)') { continue }   # refused earlier
            $full = Resolve-InputPath $v
            if ($seen.ContainsKey($full)) { continue }
            $seen[$full] = $true
            if (Test-Path -LiteralPath $full -PathType Leaf) {
                $total += [double](Get-Item -LiteralPath $full).Length
                $count++
            }
        }
    }
    $root = [System.IO.Path]::GetPathRoot($script:ProjectDir)
    $freeB = [double](New-Object System.IO.DriveInfo -ArgumentList $root).AvailableFreeSpace
    $factor = Get-Requirement 'disk_factor' 5.6
    $resultsFactor = Get-Requirement 'disk_results_factor' 1.0
    $hardMinGB = Get-Requirement 'disk_hard_min_gb' 20
    $freeGB = $freeB / 1GB
    $rootName = $root.TrimEnd('\')
    $onVolume = Get-ConfigBool 'work_volume' $true
    # With the work folder in Docker's volume, this drive holds only the results; the work files go to Docker's disk.
    if ($onVolume) {
        $needB = $total * $resultsFactor
        $workB = $total * ($factor - $resultsFactor)
    } else {
        $needB = $total * $factor
        $workB = 0
    }
    if ($count -eq 0) {
        Warn "none of the input files named in the samplesheet were found. Check the fastq_1 / fastq_2 paths in $(Split-Path $ss -Leaf)."
    }
    if ($onVolume) {
        Say ("  Disk: {0} free on {1}. Input files: {2}; the analysis may need about {3} on this drive (results) and {4} in Docker's disk (work files)." -f (Format-Bytes $freeB), $rootName, (Format-Bytes $total), (Format-Bytes $needB), (Format-Bytes $workB))
    } else {
        Say ("  Disk: {0} free on {1}. Input files: {2}; the analysis may need about {3}." -f (Format-Bytes $freeB), $rootName, (Format-Bytes $total), (Format-Bytes $needB))
    }
    if ($freeGB -lt $hardMinGB) {
        Fail ("only $(Format-Bytes $freeB) is free on $rootName. The analysis needs at least $hardMinGB GB free to start. " +
              "Free up space (or move the project to a drive with more room), then run this again.")
    }
    if ($freeB -lt $needB) {
        if ($onVolume) {
            Warn ("$rootName has $(Format-Bytes $freeB) free, but the results may need about $(Format-Bytes $needB) " +
                  "on this drive. The work files are in the Docker work volume (see Docker's disk below), so this drive " +
                  "mainly needs room for the results. If the disk fills up the run stops and can be resumed after you make room.")
        } else {
            Warn ("$rootName has $(Format-Bytes $freeB) free, but the analysis may need about $(Format-Bytes $needB). " +
                  "If the disk fills up the run stops and can be resumed after you make room. To free space: " +
                  "run_atacseq clean deletes metadata\work (the working copies; the next run rebuilds them). " +
                  "keep_work: false in config.yaml deletes them by itself after a successful run.")
        }
    }
    Check-DockerDisk $workB $onVolume
}

# Everything the run reads from the project, checked before anything is started.
function Check-ProjectInputs {
    Info "Checking the input files named in $($script:ProjectDir)\$($script:ConfigName) and samplesheet.csv"
    Check-OutputDirs
    $letters = @(Check-InputDrives)
    Check-Samplesheet
    foreach ($l in $letters) {
        Say "  ${l}:\ is mounted read-only inside the analysis (as /mnt/$($l.ToLowerInvariant())); nothing is copied from it."
    }
    if ($letters.Count -eq 0) { Say '  All inputs are inside the project folder.' }
    Check-ConfigFiles
    Check-WorkVolume
    Check-Disk
    return $letters
}

# Where Nextflow's work folder goes (config.yaml work_volume). Fails on a publish_mode
# that cannot work with the volume. Does not call Docker (it also runs in -WhatIfWsl).
function Check-WorkVolume {
    $on = Get-ConfigBool 'work_volume' $true
    $mode = Get-YamlValue $script:ConfigFile 'publish_mode'
    if (Test-Unset $mode) { $mode = '' }
    $problem = Get-PublishModeProblem $mode $on
    if ($problem) { Fail $problem }
    if ($on) {
        Say "  Work folder: Docker volume $(Get-WorkVolumeName $script:ProjectDir) (on Docker's disk, fast; resume survives restarts)"
    } else {
        Say '  Work folder: metadata\work in the project folder (work_volume: false)'
    }
}

# ---------------------------------------------------------------------------
# Project identity / container naming
# ---------------------------------------------------------------------------
function Get-ContainerName([string]$dir) {
    $leaf = (Split-Path $dir -Leaf) -replace '[^A-Za-z0-9_.-]', '_'
    $md5 = [System.Security.Cryptography.MD5]::Create()
    $h = $md5.ComputeHash([System.Text.Encoding]::UTF8.GetBytes($dir.ToLowerInvariant()))
    $short = (($h[0..3] | ForEach-Object { $_.ToString('x2') }) -join '')
    return "${PipelineName}_${leaf}_${short}"
}

function Test-Running([string]$name) {
    $out = & docker ps --filter "name=^/$name$" --format '{{.Names}}' 2>$null
    return ("$out".Trim() -eq $name)
}

# The project's work volume exists (false when Docker is down).
function Test-WorkVolume([string]$vol) {
    if (-not (Test-DockerUp)) { return $false }
    return ((Quiet { docker volume inspect $vol }) -eq 0)
}

# Bytes used by the work volume, or $null when unknown (Docker down, no image, no volume).
function Get-WorkVolumeBytes([string]$vol) {
    if (-not (Test-WorkVolume $vol)) { return $null }
    if ((Quiet { docker image inspect $Image }) -ne 0) { return $null }
    $out = @(& docker run --rm --memory 256m --mount "type=volume,source=$vol,target=/w,readonly" $Image du -sk /w 2>$null)
    foreach ($l in $out) {
        if ("$l" -match '^\s*(\d+)\s') { return ([double]$Matches[1]) * 1024 }
    }
    return $null
}

# Status lines for the work folder (used by status).
function Show-WorkFolderStatus {
    if (-not (Get-ConfigBool 'work_volume' $true)) {
        Say '  Work folder: metadata\work (work_volume: false)'
        return
    }
    if (-not (Test-DockerUp)) { Say '  Work volume: not checked (Docker is not running)'; return }
    $vol = Get-WorkVolumeName $script:ProjectDir
    if (-not (Test-WorkVolume $vol)) { Say '  Work volume: none yet'; return }
    $bytes = Get-WorkVolumeBytes $vol
    if ($null -eq $bytes) { Say "  Work volume: $vol, size not known" }
    else { Say "  Work volume: $vol, size $(Format-Bytes $bytes)" }
}

# Keep Windows from going to sleep for as long as the container runs: a hidden
# PowerShell asks for "system required" every minute and exits with the run.
# 'SSD', 'HDD', or '' (unknown: network share, USB stick, no permission...).
function Get-DriveMediaType([string]$path) {
    try {
        $letter = ([System.IO.Path]::GetPathRoot($path)).Substring(0, 1)
        $part = Get-Partition -DriveLetter $letter -ErrorAction Stop
        $disk = Get-PhysicalDisk -ErrorAction Stop | Where-Object { "$($_.DeviceId)" -eq "$($part.DiskNumber)" }
        return "$($disk.MediaType)"
    } catch { return '' }
}

# The drive letters of the FASTQ files named in the samplesheet (fastq_1 / fastq_2), e.g. 'D'. Silent on failure.
function Get-FastqDriveLetters {
    $letters = @()
    try {
        foreach ($row in @(Import-Csv -LiteralPath (Get-SamplesheetFile))) {
            foreach ($col in @('fastq_1', 'fastq_2')) {
                if (-not ($row.PSObject.Properties.Name -contains $col)) { continue }
                $v = "$($row.$col)".Trim()
                if (Test-Unset $v) { continue }
                $full = Resolve-InputPath $v
                if ($full -match '^([A-Za-z]):\\') { $letters += $matches[1].ToUpperInvariant() }
            }
        }
    } catch { }
    return @($letters | Sort-Object -Unique)
}

# ATAC_HDD for the container: '1' when one of the media types is 'HDD', else '0' (unknown '' counts as 0).
function Get-AtacHddValue($media) {
    if (@($media) -contains 'HDD') { return '1' }
    return '0'
}

# ATAC_HDD for a run: the project folder's drive, the FASTQ drives, and (with the work volume) Docker's disk.
function Get-RunHddValue([bool]$onVolume) {
    $media = @(Get-DriveMediaType $script:ProjectDir)
    foreach ($l in @(Get-FastqDriveLetters)) { $media += @(Get-DriveMediaType "$($l):\") }
    if ($onVolume) {
        $r = Get-DockerDiskRoot
        if ($r) { $media += @(Get-DriveMediaType $r) }
    }
    return (Get-AtacHddValue $media)
}

# ---------------------------------------------------------------------------
# Resources: one budget per run, enforced as the container's hard limits.
#
#   memory = 75% of the computer's RAM, but always leaving 3 GB to Windows
#            (a plain percentage starves small laptops: 80% of 8 GB leaves 1.6 GB)
#   CPUs   = all logical cores but 1 (but 2 from 8 cores up), so the desktop
#            stays responsive; heavy CPU use only slows a computer, it can't harm it
#   ... and never more than Docker Desktop itself may use (its memory limit).
#   config.yaml max_cpus / max_memory_gb replace the default (up to Docker's limit).
#
# docker run --memory makes the budget a hard ceiling (a step that outgrows it
# is killed inside the container and retried, Windows never starts swapping);
# --cpuset-cpus makes tools that use "all cores" (cutadapt --cores 0) see
# only the budget's cores.
# ---------------------------------------------------------------------------
function Get-Requirement([string]$key, [double]$default) {
    $f = Join-Path $PipelineDir 'docker\requirements.conf'
    if (Test-Path $f) {
        foreach ($line in [System.IO.File]::ReadAllLines($f)) {
            if ($line -match "^\s*$([regex]::Escape($key))\s*=\s*(\d+(\.\d+)?)") { return [double]::Parse($Matches[1], [System.Globalization.CultureInfo]::InvariantCulture) }
        }
    }
    return $default
}

function Get-HostResources {
    $r = @{ memMB = 0; cpus = 1; freeMB = 0 }
    try {
        $cs = Get-CimInstance Win32_ComputerSystem
        $r.memMB = [int][math]::Floor($cs.TotalPhysicalMemory / 1MB)
        $r.cpus  = [int]$cs.NumberOfLogicalProcessors
        $r.freeMB = [int][math]::Floor((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1024)
    } catch { }
    return $r
}

function Get-DockerResources {
    $r = @{ memMB = 0; cpus = 0 }
    $n = "$(& docker info --format '{{.NCPU}}' 2>$null)".Trim()
    $m = "$(& docker info --format '{{.MemTotal}}' 2>$null)".Trim()
    if ($n -match '^\d+$') { $r.cpus = [int]$n }
    if ($m -match '^\d+$') { $r.memMB = [int][math]::Floor([double]$m / 1MB) }
    return $r
}

# The default policy (no config.yaml override), from the computer alone.
function Get-PolicyBudget($h) {
    $mem = [math]::Min([math]::Floor($h.memMB * 0.75), $h.memMB - 3072)
    $spare = 1; if ($h.cpus -ge 8) { $spare = 2 }
    $cpus = [math]::Max(1, $h.cpus - $spare)
    return @{ memMB = [int][math]::Max(0, $mem); cpus = [int]$cpus }
}

# Returns @{ cpus; memMB; reason } for this run, or stops with an explanation.
function Get-RunBudget {
    $h = Get-HostResources
    $d = Get-DockerResources
    $p = Get-PolicyBudget $h
    $dockMem = $d.memMB - 512           # the Docker VM's own kernel and services
    if ($dockMem -le 0) { $dockMem = $p.memMB }
    $dockCpu = $d.cpus; if ($dockCpu -le 0) { $dockCpu = $p.cpus }

    $wantMem = Get-YamlValue $script:ConfigFile 'max_memory_gb'
    $wantCpu = Get-YamlValue $script:ConfigFile 'max_cpus'
    $why = @()
    if ("$wantMem" -match '^\d+(\.\d+)?$') {
        $wantGB = [double]::Parse("$wantMem", [System.Globalization.CultureInfo]::InvariantCulture)
        $mem = [int][math]::Min($wantGB * 1024, $dockMem)
        $why += "max_memory_gb: $wantMem in $($script:ConfigName)"
        if ($wantGB * 1024 -gt $dockMem) { Warn "max_memory_gb: $wantMem is more than Docker Desktop may use; using $(Format-GB $dockMem)." }
    } else {
        $mem = [int][math]::Min($p.memMB, $dockMem)
        $useful = [int][math]::Round((Get-Requirement 'max_useful_memory_gb' 0) * 1024)
        if ($useful -gt 0 -and $mem -gt $useful) {
            $mem = $useful; $why += "the most this pipeline can use (docker\requirements.conf)"
        }
        elseif ($p.memMB -le $dockMem) { $why += "75% of this computer's $(Format-GB $h.memMB), leaving at least 3 GB to Windows" }
        else { $why += "all that Docker Desktop may use: $(Format-GB $d.memMB) of this computer's $(Format-GB $h.memMB)" }
    }
    if ("$wantCpu" -match '^\d+$' -and [int]$wantCpu -ge 1) {
        $cpus = [int][math]::Min([int]$wantCpu, $dockCpu)
        $why += "max_cpus: $wantCpu"
    } else {
        $cpus = [int][math]::Min($p.cpus, $dockCpu)
    }

    $minMem = [int][math]::Round((Get-Requirement 'min_memory_gb' 2.5) * 1024)
    $minCpu = [int](Get-Requirement 'min_cpus' 1)
    if ($mem -lt $minMem) {
        if ("$wantMem" -match '^\d') {
            Fail "max_memory_gb: $wantMem in $($script:ConfigName) is below the about $(Format-GB $minMem) this pipeline needs. Raise it or remove the line."
        } elseif ($p.memMB -ge $minMem) {
            Fail ("Docker Desktop may only use $(Format-GB $d.memMB), but this pipeline needs about $(Format-GB $minMem) " +
                  "(plus Docker's own share). Raise Docker Desktop's memory limit (Settings > Resources > Memory), " +
                  "then run this again.")
        } else {
            Fail ("This computer has $(Format-GB $h.memMB) of memory. Keeping 3 GB for Windows leaves less than the about " +
                  "$(Format-GB $minMem) this pipeline needs. Use a computer with at least $([math]::Ceiling(($minMem + 3072) / 1024)) GB, " +
                  "or set max_memory_gb in config.yaml to use more (Windows may then become slow while it runs).")
        }
    }
    if ($cpus -lt $minCpu) { Warn "only $cpus CPU core(s) for the analysis: it will work, but slowly." }
    if ($h.freeMB -gt 0 -and $h.freeMB -lt $mem) {
        Warn ("only $(Format-GB $h.freeMB) of memory is free right now (other programs use the rest), and the analysis may use up to " +
              "$(Format-GB $mem). Closing large programs (e.g. a browser with many tabs) keeps Windows responsive.")
    }
    return @{ cpus = $cpus; memMB = $mem; reason = ($why -join '; ') }
}

# ---------------------------------------------------------------------------
# Docker's memory limit for an analysis: the memory= line of [wsl2] in
# %USERPROFILE%\.wslconfig (Docker Desktop's VM). The functions from
# Split-TextLines to Remove-WslMemoryText work on text only (no files, no
# Docker), so they can be tested on sample text (see selftest).
# ---------------------------------------------------------------------------
function New-TextLine([string]$t, [string]$e) { return @{ T = $t; E = $e } }

# Lines with their own line ending ("", "\n" or "\r\n"), so untouched lines
# are written back byte for byte.
function Split-TextLines([string]$text) {
    $out = New-Object System.Collections.ArrayList
    $start = 0
    for ($i = 0; $i -lt $text.Length; $i++) {
        if ($text[$i] -eq [char]10) {
            $end = $i; $eol = "`n"
            if ($i -gt $start -and $text[$i - 1] -eq [char]13) { $end = $i - 1; $eol = "`r`n" }
            [void]$out.Add((New-TextLine $text.Substring($start, $end - $start) $eol))
            $start = $i + 1
        }
    }
    if ($start -lt $text.Length) { [void]$out.Add((New-TextLine $text.Substring($start) '')) }
    return $out.ToArray()
}

function Join-TextLines($lines) {
    $sb = New-Object System.Text.StringBuilder
    foreach ($l in $lines) { [void]$sb.Append([string]$l.T); [void]$sb.Append([string]$l.E) }
    return $sb.ToString()
}

function Get-NewlineStyle([string]$text) {
    if ($text.Contains("`r`n")) { return "`r`n" }
    if ($text.Contains("`n")) { return "`n" }
    return "`r`n"
}

function Get-HeaderName([string]$t) {
    if ($t -match '^\s*\[([^\]]+)\]\s*$') { return $Matches[1].Trim() }
    return $null
}

# memory= line: Prefix ("memory = "), Value ("6GB"), Rest (spaces and any comment).
function Get-MemoryParts([string]$t) {
    if ($t -match '^(\s*memory\s*=\s*)([^#;\s]*)(.*)$') {
        return @{ Prefix = $Matches[1]; Value = $Matches[2]; Rest = $Matches[3] }
    }
    return $null
}

# The old value recorded by an earlier run: "6GB", "none" or "none (section created)".
function Get-MarkerValue([string]$t) {
    if ($t -match '^\s*#\s*run_atacseq: previous memory=(.*?)\s*$') { return $Matches[1] }
    return $null
}

function Find-MarkerIndex($lines) {
    for ($i = 0; $i -lt $lines.Count; $i++) {
        if ($null -ne (Get-MarkerValue $lines[$i].T)) { return $i }
    }
    return -1
}

# First [wsl2] section: Header = its line index (-1 if none), End = index where it ends.
function Find-Wsl2Section($lines) {
    $hdr = -1; $end = $lines.Count
    for ($i = 0; $i -lt $lines.Count; $i++) {
        $n = Get-HeaderName $lines[$i].T
        if ($null -ne $n) {
            if ($hdr -ge 0) { $end = $i; break }
            if ($n -ieq 'wsl2') { $hdr = $i }
        }
    }
    return @{ Header = $hdr; End = $end }
}

function Find-MemoryLine($lines, $sec) {
    if ($sec.Header -lt 0) { return -1 }
    for ($i = $sec.Header + 1; $i -lt $sec.End; $i++) {
        if ($null -ne (Get-MemoryParts $lines[$i].T)) { return $i }
    }
    return -1
}

function Get-WslMemoryValue([string]$text) {
    $lines = @(Split-TextLines $text)
    $sec = Find-Wsl2Section $lines
    $mi = Find-MemoryLine $lines $sec
    if ($mi -lt 0) { return $null }
    return (Get-MemoryParts $lines[$mi].T).Value
}

# "6GB", "8 GB", "512MB" -> GB as a number; $null when it cannot be read.
function Convert-WslValueToGB([string]$v) {
    if ($v -match '^\s*(\d+(?:\.\d+)?)\s*(GB|G|MB|M)\s*$') {
        $n = [double]::Parse($Matches[1], [System.Globalization.CultureInfo]::InvariantCulture)
        if ($Matches[2] -match '^(MB|M)$') { $n = $n / 1024.0 }
        return $n
    }
    return $null
}

# The analysis setting by total RAM: under 12 GB -> 6 GB, 12 to 24 GB -> 8 GB, 24 GB or more -> 16 GB.
function Get-WslTierGB([double]$ramGB) {
    $ramGB = [math]::Round($ramGB)   # Windows reports a little less than the nominal size (a 12 GB laptop shows 11.7)
    if ($ramGB -lt 12) { return 6 }
    if ($ramGB -lt 24) { return 8 }
    return 16
}

# Docker's limit now: the memory= line; else Docker's own value (docker info) when
# it runs; else Docker's default (half of the computer's memory).
function Get-WslCurrentGB([string]$text, [double]$ramGB, [int]$dockerMB) {
    $v = Get-WslMemoryValue $text
    if ($null -ne $v) {
        $g = Convert-WslValueToGB $v
        if ($null -ne $g) { return [pscustomobject]@{ GB = [int][math]::Round($g); Source = 'memory= in .wslconfig' } }
    }
    if ($dockerMB -gt 0) { return [pscustomobject]@{ GB = [int][math]::Round($dockerMB / 1024.0); Source = "Docker's current limit" } }
    return [pscustomobject]@{ GB = [int][math]::Round($ramGB / 2.0); Source = "Docker's default, half of the computer's memory" }
}

# Adds memory=<gb>GB to [wsl2] (the first one; a new section is put at the top
# of the file). The line gets a marker above it with the old value.
function Add-WslMemoryText([string]$text, [int]$gb) {
    $nl = Get-NewlineStyle $text
    $lines = @(Split-TextLines $text)
    $val = "${gb}GB"
    $sec = Find-Wsl2Section $lines
    $out = New-Object System.Collections.ArrayList
    if ($sec.Header -ge 0) {
        $mi = Find-MemoryLine $lines $sec
        if ($mi -ge 0) {
            $p = Get-MemoryParts $lines[$mi].T
            for ($i = 0; $i -lt $lines.Count; $i++) {
                if ($i -eq $mi) {
                    [void]$out.Add((New-TextLine "# run_atacseq: previous memory=$($p.Value)" $nl))
                    [void]$out.Add((New-TextLine ($p.Prefix + $val + $p.Rest) $lines[$i].E))
                } else {
                    [void]$out.Add($lines[$i])
                }
            }
        } else {
            for ($i = 0; $i -lt $lines.Count; $i++) {
                [void]$out.Add($lines[$i])
                if ($i -eq $sec.Header) {
                    $memEol = $nl
                    if ($lines[$i].E -eq '') { $lines[$i].E = $nl; $memEol = '' }   # keep "no newline at the end of the file"
                    [void]$out.Add((New-TextLine '# run_atacseq: previous memory=none' $nl))
                    [void]$out.Add((New-TextLine "memory=$val" $memEol))
                }
            }
        }
    } else {
        [void]$out.Add((New-TextLine '[wsl2]' $nl))
        [void]$out.Add((New-TextLine '# run_atacseq: previous memory=none (section created)' $nl))
        [void]$out.Add((New-TextLine "memory=$val" $nl))
        for ($i = 0; $i -lt $lines.Count; $i++) { [void]$out.Add($lines[$i]) }
    }
    return (Join-TextLines $out.ToArray())
}

# Undoes Add-WslMemoryText. Found = $false when there is no marker (text unchanged).
function Remove-WslMemoryText([string]$text) {
    $lines = @(Split-TextLines $text)
    $mi = Find-MarkerIndex $lines
    if ($mi -lt 0) { return @{ Text = $text; Found = $false; Old = $null; Created = $false } }
    $old = Get-MarkerValue $lines[$mi].T
    $created = ($old -eq 'none (section created)')
    $none = ($old -eq 'none')
    $drop = @{}
    $drop[$mi] = $true
    $fixIdx = -1
    if (($mi + 1) -lt $lines.Count -and $null -ne (Get-MemoryParts $lines[$mi + 1].T)) {
        if ($created -or $none) { $drop[$mi + 1] = $true } else { $fixIdx = $mi + 1 }
    }
    if ($created -and $mi -gt 0 -and ((Get-HeaderName $lines[$mi - 1].T) -ieq 'wsl2')) { $drop[$mi - 1] = $true }
    $out = New-Object System.Collections.ArrayList
    for ($i = 0; $i -lt $lines.Count; $i++) {
        if ($drop.ContainsKey($i)) { continue }
        if ($i -eq $fixIdx) {
            $p = Get-MemoryParts $lines[$i].T
            [void]$out.Add((New-TextLine ($p.Prefix + $old + $p.Rest) $lines[$i].E))
        } else {
            [void]$out.Add($lines[$i])
        }
    }
    return @{ Text = (Join-TextLines $out.ToArray()); Found = $true; Old = $old; Created = $created }
}

# The decision for one run, as text: Action = change | none | marked | busy | optout.
function Get-WslPlan([string]$text, [double]$ramGB, [int]$dockerMB, [bool]$busy, [bool]$optOut) {
    $tier = Get-WslTierGB $ramGB
    $cur = Get-WslCurrentGB $text $ramGB $dockerMB
    $plan = [pscustomobject]@{ Tier = $tier; Current = $cur; Action = 'none'; NewText = $text }
    if ($optOut) { $plan.Action = 'optout'; return $plan }
    if ((Find-MarkerIndex @(Split-TextLines $text)) -ge 0) { $plan.Action = 'marked'; return $plan }
    if ($cur.GB -eq $tier) { return $plan }
    if ($busy) { $plan.Action = 'busy'; return $plan }
    $plan.Action = 'change'
    $plan.NewText = Add-WslMemoryText $text $tier
    return $plan
}

function Show-WslDiff([string]$old, [string]$new) {
    $oldLines = @(Split-TextLines $old)
    $newLines = @(Split-TextLines $new)
    $oldSet = @{}; foreach ($l in $oldLines) { $oldSet[$l.T] = $true }
    $newSet = @{}; foreach ($l in $newLines) { $newSet[$l.T] = $true }
    Say '  Lines removed:'
    $any = $false
    foreach ($l in $oldLines) { if (-not $newSet.ContainsKey($l.T)) { Say "    - $($l.T)"; $any = $true } }
    if (-not $any) { Say '    (none)' }
    Say '  The file after the change (added lines marked with +):'
    if ($newLines.Count -eq 0) { Say '    (empty)' }
    foreach ($l in $newLines) {
        if ($oldSet.ContainsKey($l.T)) { Say "      $($l.T)" } else { Say "    + $($l.T)" }
    }
}

function Show-WslPlanReport($plan, [double]$ramGB, [string]$oldText, [bool]$dockerUp) {
    Say "Docker's memory limit for this analysis:"
    Say ("  Memory of this computer : {0:0.#} GB" -f $ramGB)
    Say "  Setting for analyses    : $($plan.Tier) GB"
    Say "  Now                     : $($plan.Current.GB) GB ($($plan.Current.Source))"
    switch ($plan.Action) {
        'none'   { Say '  Decision                : no change needed' }
        'marked' { Say '  Decision                : already set by an earlier analysis; it is put back when that one finishes' }
        'optout' { Say '  Decision                : not changed (--no-wsl-resize or wsl_resize: false in config.yaml)' }
        'busy'   { Say '  Decision                : not changed: other containers are running' }
        'change' {
            Say "  Decision                : set it to $($plan.Tier) GB for this analysis, then put $($plan.Current.GB) GB back"
            Say '  Would run               : wsl --shutdown, then restart Docker Desktop'
            Say '  Afterwards (nothing else running): write the old value back, wsl --shutdown, restart Docker Desktop'
            Show-WslDiff $oldText $plan.NewText
        }
    }
    if (-not $dockerUp) { Say '  (Docker is not running, so its current limit is taken from the file or its default. Nothing is started.)' }
}

# Writes the file (a copy of the user's file is kept once, as .wslconfig.run_atacseq.bak).
function Write-WslFile([string]$newText) {
    if ((Test-Path -LiteralPath $script:WslFile) -and -not (Test-Path -LiteralPath $script:WslBackup)) {
        Copy-Item -LiteralPath $script:WslFile -Destination $script:WslBackup
    }
    Write-TextRaw $script:WslFile $newText
}

function Restart-DockerDesktop {
    Info 'Restarting Docker Desktop to apply the memory limit (about a minute)...'
    Get-Process -Name 'Docker Desktop', 'com.docker.backend', 'com.docker.build' -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
    Start-Sleep -Seconds 3
    & wsl.exe --shutdown
    Start-Sleep -Seconds 5
    Ensure-Docker
}

# Called by run before the analysis starts (after the input checks).
function Invoke-WslStep {
    $optOut = $script:NoWslResize -or (-not (Get-ConfigBool 'wsl_resize' $true))
    $ramGB = (Get-HostResources).memMB / 1024.0
    $text = Read-TextRaw $script:WslFile
    if ($null -eq $text) { $text = '' }
    $dockerUp = Test-DockerUp
    $dockerMB = 0; $busy = 0
    if ($dockerUp) { $dockerMB = (Get-DockerResources).memMB; $busy = Get-RunningCount }
    $plan = Get-WslPlan $text $ramGB $dockerMB ($busy -gt 0) $optOut

    if ($script:WhatIfWsl) { Show-WslPlanReport $plan $ramGB $text $dockerUp; return }

    switch ($plan.Action) {
        'change' {
            Say ''
            Say "To run this analysis safely, Docker's memory limit needs to be $($plan.Tier) GB (it is now $($plan.Current.GB) GB)."
            Say "I'll set it, restart Docker (about a minute), and put the old value back when the analysis is done."
            if (-not (Confirm "Set Docker's memory limit to $($plan.Tier) GB for this analysis?")) {
                Say "Left at $($plan.Current.GB) GB: the analysis runs with what Docker has now."
                return
            }
            Write-WslFile $plan.NewText
            $script:WslRestorePending = $true     # from here on the old value is put back, whatever happens
            Restart-DockerDesktop
            Say "Docker's memory limit is now $($plan.Tier) GB. It goes back to its previous setting when the analysis is done."
        }
        'busy' {
            Say "Docker is running other containers right now, so its memory limit stays at $($plan.Current.GB) GB for this analysis."
            Say '  (Stop them and run again to let the analysis use the setting for analyses.)'
        }
        'marked' {
            Say "Docker's memory limit is still set for an earlier analysis; it is put back when that one finishes."
        }
    }
}

# Puts the old value back when a marker exists and nothing else is running.
function Restore-WslIfIdle([switch]$Ask) {
    if ($script:WhatIfWsl) { return }
    $text = Read-TextRaw $script:WslFile
    if ($null -eq $text) { return }
    $lines = @(Split-TextLines $text)
    $mi = Find-MarkerIndex $lines
    if ($mi -lt 0) { return }
    $old = Get-MarkerValue $lines[$mi].T
    if ($old -eq 'none' -or $old -eq 'none (section created)') { $label = "Docker's default (no memory line)" } else { $label = "$old" }
    if (-not (Test-DockerUp)) {
        Say "Docker's memory limit was set for an analysis that was closed early. Docker is not running, so it cannot be put back now;"
        Say "it is put back the next time a run_atacseq command runs while Docker is running."
        return
    }
    if ((Get-RunningCount) -gt 0) {
        Say "Docker's memory limit stays at its analysis setting for now, because other containers are running."
        Say "It is put back by the next run_atacseq command that finds nothing running."
        return
    }
    if ($Ask) {
        if (-not (Confirm "Docker's memory limit was changed for an analysis that was closed early. Put it back to $label now? (Docker restarts, about a minute)")) {
            Say 'Left as it is.'
            return
        }
    }
    $r = Remove-WslMemoryText $text
    if ($r.Created -and $r.Text.Length -eq 0) {
        Remove-Item -LiteralPath $script:WslFile -Force
    } else {
        Write-TextRaw $script:WslFile $r.Text
    }
    if (Test-Path -LiteralPath $script:WslBackup) { Remove-Item -LiteralPath $script:WslBackup -Force }
    Restart-DockerDesktop
    Say "Docker's memory limit is back to $label."
}

# ---------------------------------------------------------------------------
# Container arguments and log following
# ---------------------------------------------------------------------------
# Keep Windows awake and guard its memory while a container (or an image
# download) runs: Docker's Linux VM keeps the files it has read as cache, and
# hands that memory back slowly; when Windows' free memory drops below
# max(4 GB, 20% of RAM) the helper tells the VM to drop that cache.
function Start-Guard([string]$container = '', [int]$watchPid = 0) {
    $totalMB = 0
    try { $totalMB = [int]((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1MB) } catch { }
    $minFreeMB = [math]::Max(4096, [int]($totalMB / 5))
    $alive = "[bool](Get-Process -Id $watchPid -ErrorAction SilentlyContinue)"
    if ($container) { $alive = "[bool](& docker ps -q --filter `"name=^/$container`$`" 2>`$null)" }
    $code = @"
Add-Type -Namespace AtacSeqPower -Name Power -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint flags);'
`$i = 0
while ($alive) {
    if (`$i % 6 -eq 0) { [AtacSeqPower.Power]::SetThreadExecutionState([uint32]2147483649) | Out-Null }
    `$freeMB = 1e9
    try { `$freeMB = (Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1024 } catch { }
    if (`$freeMB -lt $minFreeMB) { & wsl.exe -d docker-desktop -u root sh -c 'echo 1 > /proc/sys/vm/drop_caches' 2>`$null | Out-Null }
    `$i++
    Start-Sleep -Seconds 10
}
"@
    $encoded = [Convert]::ToBase64String([System.Text.Encoding]::Unicode.GetBytes($code))
    try {
        Start-Process powershell.exe -WindowStyle Hidden -ArgumentList '-NoProfile', '-ExecutionPolicy', 'Bypass', '-EncodedCommand', $encoded | Out-Null
    } catch { }
}

function Follow-Run([string]$name) {
    & docker logs -f $name 2>&1 | ForEach-Object { "$_" }
    # The container is gone once the run ends; its outcome is in the project.
    Start-Sleep -Seconds 1
    $statusFile = Join-Path $script:ProjectDir 'metadata\nf_logs\last_run.status'
    $st = ''
    if (Test-Path $statusFile) { $st = ([System.IO.File]::ReadAllText($statusFile)).Trim() }
    Say ''
    if ($st -eq 'succeeded') {
        Info 'Pipeline finished successfully.'
        Say  "  Results : $($script:ProjectDir)\results\"
        # keep_work: false: the container has emptied the work volume; remove the volume too.
        # (With --detach nothing runs here: an empty volume stays until the next clean.)
        if ((-not (Get-ConfigBool 'keep_work' $true)) -and (Get-ConfigBool 'work_volume' $true)) {
            $vol = Get-WorkVolumeName $script:ProjectDir
            if (Test-WorkVolume $vol) {
                if ((Quiet { docker volume rm $vol }) -eq 0) { Say "  The Docker work volume $vol was deleted (keep_work: false)." }
            }
        }
        Done 0
    } elseif ($st -eq 'running' -and (Test-Running $name)) {
        Say 'Stopped watching; the analysis is still running.'
        Done 0
    } else {
        Write-Host "Pipeline did not finish: $st" -ForegroundColor Red
        Say '  See:  run_atacseq status   (it names the failed step and where its messages were saved)'
        Say '  Fix the problem, then run again: finished steps are reused, not redone.'
        Done 1
    }
}

# Drives named by script arguments (e.g. utils\X.R D:\data\y.csv) are mounted too.
function Add-ArgDrives([string[]]$letters, $argList) {
    $all = @($letters)
    foreach ($a in @($argList)) {
        $s = [string]$a
        if ($s -match '^(\\\\|//)') { Fail "A script argument names a network path ($s). Map the share to a drive letter first." }
        if ($s -match '^([A-Za-z]):[\\/]') {
            $l = $Matches[1].ToUpperInvariant()
            if (-not (Test-Path -LiteralPath "$($l):\")) { Fail "The script argument $s needs the drive $($l):, which is not connected." }
            $all += $l
        }
    }
    return @($all | Sort-Object -Unique)
}

# Check the project is not inside this pipeline's own repository folder.
function Check-NotInsidePipeline([string]$dir) {
    $repoRoot = (Split-Path $PipelineDir -Parent).TrimEnd('\')
    if ($dir.Equals($repoRoot, [StringComparison]::OrdinalIgnoreCase) -or $dir.StartsWith($repoRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
        Fail "Don't put a project inside the pipeline's own folder ($repoRoot). Choose another folder, e.g. $HOME\Documents\my_atacseq"
    }
}

# Path relative to the project folder with forward slashes, or $null if outside it.
function Get-ProjectRelative([string]$full) {
    $root = $script:ProjectDir.TrimEnd('\') + '\'
    if ($full.Equals($script:ProjectDir.TrimEnd('\'), [StringComparison]::OrdinalIgnoreCase)) { return '.' }
    if ($full.StartsWith($root, [StringComparison]::OrdinalIgnoreCase)) {
        return ($full.Substring($root.Length) -replace '\\', '/')
    }
    return $null
}

function Get-ProjectDirFromArgs {
    $script:ProjectDir = (Get-Location).ProviderPath
    $script:ConfigName = 'config.yaml'
    if ($positional.Count -gt 0) {
        $c = $positional[0]
        if (Test-Path -LiteralPath $c -PathType Container) {
            $script:ProjectDir = (Resolve-Path -LiteralPath $c).ProviderPath
        } elseif (Test-Path -LiteralPath $c -PathType Leaf) {
            $item = Get-Item -LiteralPath $c
            $script:ProjectDir = $item.DirectoryName
            $script:ConfigName = $item.Name
        } else {
            Fail "$c not found (expected a config .yaml file or a project folder)."
        }
    }
    $script:ProjectDir = $script:ProjectDir.TrimEnd('\')
    $script:ConfigFile = Join-Path $script:ProjectDir $script:ConfigName
}

# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------
function Show-Help {
    Say @"
ATAC-seq pipeline (Docker version for Windows)

  run_atacseq setup            one-time: Docker Desktop + the pipeline's software image
  run_atacseq init [FOLDER]    set up a new project folder (default: current folder)
  run_atacseq                  run the pipeline for the project in the current folder
  run_atacseq run CONFIG.yaml  ... with another config file
  run_atacseq status           progress of the current/last run (and why it failed)
  run_atacseq logs             follow the live log
  run_atacseq stop             stop the running analysis (re-running resumes it)
  run_atacseq clean            delete the work folder, resume cache and Docker work volume now (results\ is kept)
  run_atacseq script utils\X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment

Options: --detach  --yes  --pause  --no-wsl-resize (run: keep Docker's memory limit as it is)
script: foreground only (no --detach); arguments with double quotes are not supported on Windows.
Pipeline folder: $PipelineDir
"@
}

function Cmd-Setup {
    Info 'ATAC-seq setup (Docker version for Windows)'
    Ensure-Docker -AllowInstall
    Repair-LineEndings
    Pull-Image

    # What one analysis may use on this computer (the run sets Docker's limit to fit it).
    $h = Get-HostResources
    $d = Get-DockerResources
    $p = Get-PolicyBudget $h
    $tierMB = (Get-WslTierGB ($h.memMB / 1024.0)) * 1024
    $mem = [math]::Min($p.memMB, $tierMB - 512)
    $useful = [int][math]::Round((Get-Requirement 'max_useful_memory_gb' 0) * 1024)
    if ($useful -gt 0) { $mem = [math]::Min($mem, $useful) }
    Say ''
    Say "This computer: $($h.cpus) CPU cores, $(Format-GB $h.memMB) memory. Docker Desktop may use now: $($d.cpus) cores, $(Format-GB $d.memMB)."
    Say "When an analysis starts, run_atacseq sets Docker's memory limit to suit this computer (it asks first) and puts it back afterwards."
    Say "Each analysis then uses at most $($p.cpus) CPU cores and about $(Format-GB $mem) of memory (max_cpus / max_memory_gb in config.yaml change this)."
    $minGB = Get-Requirement 'min_memory_gb' 2.5
    if ($mem -lt $minGB * 1024) {
        Warn "that is less than the about $minGB GB this pipeline needs: runs will refuse to start."
    }

    # Tool versions, straight from the image.
    Say ''
    # No double quotes inside the bash command: Windows PowerShell 5.1 would strip them.
    & docker run --rm --platform linux/amd64 $Image bash -c 'echo nextflow: $(nextflow -version 2>&1 | grep -i version | head -1); echo R: $(R --version 2>&1 | head -1); echo python3: $(python3 --version 2>&1); echo samtools: $(samtools --version 2>&1 | head -1); echo bowtie2: $(bowtie2 --version 2>&1 | head -1)'

    # Put run_atacseq on the user's PATH (new terminals only).
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if (-not $userPath) { $userPath = '' }
    $parts = $userPath -split ';' | Where-Object { $_ }
    if (-not ($parts | Where-Object { $_.TrimEnd('\') -ieq $PipelineDir.TrimEnd('\') })) {
        [Environment]::SetEnvironmentVariable('Path', (($parts + $PipelineDir) -join ';'), 'User')
        Say ''
        Say "Added $PipelineDir to your PATH: in a NEW terminal window, 'run_atacseq' works from any folder."
    }
    Say ''
    Info 'Setup complete. Next: create a project folder with   run_atacseq init <folder>'
    Done 0
}

function Cmd-Init {
    $target = (Get-Location).ProviderPath
    if ($positional.Count -gt 0) { $target = $positional[0] }
    New-Item -ItemType Directory -Force $target | Out-Null
    $target = (Resolve-Path -LiteralPath $target).ProviderPath.TrimEnd('\')
    Check-NotInsidePipeline $target
    Info "Setting up a new $PipelineName project in $target ..."
    foreach ($d in @('original_data\fastq', 'metadata', 'results', 'utils')) {
        New-Item -ItemType Directory -Force (Join-Path $target $d) | Out-Null
    }
    $cfg = Join-Path $target 'config.yaml'
    if (Test-Path $cfg) { Say '  config.yaml already exists - leaving it alone.' }
    else { Copy-Item (Join-Path $PipelineDir 'examples\config.yaml') $cfg; Say '  Fetched config.yaml' }

    # CSV templates: copied into examples\ (never overwritten). Copy the ones you need to the project folder and fill them in.
    $ex = Join-Path $target 'examples'
    New-Item -ItemType Directory -Force $ex | Out-Null
    $csvCount = 0
    foreach ($f in Get-ChildItem (Join-Path $PipelineDir 'examples') -Filter '*.csv' -File) {
        $d = Join-Path $ex $f.Name
        if (-not (Test-Path $d)) { Copy-Item $f.FullName $d; $csvCount++ }
    }
    Say "  Fetched examples\ ($csvCount CSV templates, e.g. samplesheet.csv)"

    # AI agent instructions and agents, when this pipeline folder has them (agent\ holds them).
    $srcAgent = Join-Path $PipelineDir 'agent'
    foreach ($f in @('AGENTS.md', 'CLAUDE.md')) {
        $src = Join-Path $srcAgent $f
        if (-not (Test-Path -LiteralPath $src)) { continue }
        $d = Join-Path $target $f
        if (Test-Path $d) { Say "  $f already exists - leaving it alone."; continue }
        $text = [System.IO.File]::ReadAllText($src)
        Write-TextLF $d ($text.Replace('{{PIPELINE_DIR}}', $PipelineDir))
        Say "  Fetched $f"
    }
    $srcAgents = Join-Path $srcAgent '.claude\agents'
    if (Test-Path -LiteralPath $srcAgents -PathType Container) {
        $dstAgents = Join-Path $target '.claude\agents'
        $srcRoot = $srcAgents.TrimEnd('\') + '\'
        $copied = 0
        foreach ($f in Get-ChildItem -LiteralPath $srcAgents -Recurse -File -Force) {
            $rel = $f.FullName.Substring($srcRoot.Length)
            $dst = Join-Path $dstAgents $rel
            if (Test-Path -LiteralPath $dst) { continue }
            New-Item -ItemType Directory -Force (Split-Path $dst -Parent) | Out-Null
            Copy-Item -LiteralPath $f.FullName -Destination $dst
            $copied++
        }
        if ($copied -gt 0) { Say "  Fetched .claude\agents\ ($copied files for AI agents; existing files left alone)" }
    }

    # Launchers inside the project (rewritten each time: they hold the pipeline's location).
    # The terminal launcher lives in utils\ (the project root keeps only the double-click RUN_ANALYSIS.bat).
    $bat = Join-Path $PipelineDir 'run_atacseq.bat'
    $marker = 'rem Runs the ATAC-seq pipeline (Docker version) for THIS project folder.'
    $cmd = "@echo off`r`n$marker`r`n" +
           "rem Same commands as the pipeline's own launcher: (no argument) = run, status, logs, stop, clean, script utils\X.R [ARGS].`r`n" +
           "rem Generated by 'run_atacseq init'; the pipeline lives in:`r`nrem   $PipelineDir`r`n" +
           "cd /d `"%~dp0..`"`r`ncall `"$bat`" %*`r`nexit /b %ERRORLEVEL%`r`n"
    [System.IO.File]::WriteAllText((Join-Path $target 'utils\run_atacseq.bat'), $cmd, [System.Text.Encoding]::ASCII)
    $dbl = "@echo off`r`nrem Double-click to run (or resume) the ATAC-seq analysis of this folder.`r`n" +
           "cd /d `"%~dp0`"`r`ncall `"$bat`" run --pause`r`n"
    [System.IO.File]::WriteAllText((Join-Path $target 'RUN_ANALYSIS.bat'), $dbl, [System.Text.Encoding]::ASCII)
    Say '  Created utils\run_atacseq.bat (terminal) and RUN_ANALYSIS.bat (double-click)'

    # .gitattributes: LF line endings for the scripts, so a utils\ copy stays byte-identical
    # to the pipeline's default on every OS. Never overwritten.
    $ga = Join-Path $target '.gitattributes'
    if (-not (Test-Path -LiteralPath $ga)) {
        $gaLines = @(
            '# Keep scripts byte-identical to the pipeline''s copies on every OS: utils/',
            '# replacements are compared with the pipeline''s default scripts, and bash',
            '# scripts need LF line endings.',
            'utils/** text eol=lf',
            '*.bat text eol=crlf',
            '*.ps1 text eol=crlf',
            '*.sh text eol=lf',
            '*.R text eol=lf',
            '*.py text eol=lf'
        )
        Write-TextLF $ga (($gaLines -join "`n") + "`n")
        Say '  Fetched .gitattributes (LF line endings for scripts)'
    }

    if (Test-Path (Join-Path $target '.git')) {
        Say '  Already a git repository - leaving it alone.'
    } elseif (Have 'git') {
        $gi = Join-Path $target '.gitignore'
        if (-not (Test-Path $gi)) { Copy-Item (Join-Path $PipelineDir 'examples\project.gitignore') $gi }
        $null = & git -C $target init -q 2>&1
        $paths = @()
        foreach ($p in @('config.yaml', 'examples', '.gitignore', '.gitattributes', 'AGENTS.md', 'CLAUDE.md', 'RUN_ANALYSIS.bat', 'utils', '.claude')) {
            if (Test-Path -LiteralPath (Join-Path $target $p)) { $paths += $p }
        }
        $null = & git -C $target add @paths 2>&1
        $null = & git -C $target -c user.name=atacseq -c user.email=atacseq@localhost commit -q -m "Initialize $PipelineName project" 2>&1
        Say '  Initialized git repository (raw data / results / metadata excluded - see .gitignore)'
    } else {
        Say '  (git not installed - skipped version control; optional)'
    }

    Say ''
    Say 'Next steps:'
    Say '  1. Copy the input tables you need from examples\ into this folder and fill them in:'
    Say '       samplesheet.csv (one row per sample) and, if you merge peaks, peak_merges.csv (see config.yaml).'
    Say '  2. Edit config.yaml: project_name, samplesheet, and the reference files (genome_index, contamination_index, blacklist, chrom_sizes, gtf, tss_bed).'
    Say '     Write Windows paths in full, e.g. D:\data\x.fastq.gz, or relative to this folder, e.g. original_data\x.fastq.gz.'
    Say '     Network paths (\\server\share\...) are not supported: map the share to a drive letter first.'
    Say '     Your input files stay where they are: they are never moved or copied.'
    Say '  3. Run:  double-click RUN_ANALYSIS.bat   (or type  .\utils\run_atacseq.bat  in this folder)'
    Done 0
}

function Cmd-Run {
    Get-ProjectDirFromArgs
    if (-not (Test-Path -LiteralPath $script:ConfigFile)) {
        Fail "$($script:ConfigName) not found in $($script:ProjectDir). Set the project up first:  run_atacseq init"
    }
    Check-NotInsidePipeline $script:ProjectDir

    if ($script:WhatIfWsl) {
        # Maintainer dry run: reads the project and the Docker memory state, changes nothing, starts nothing.
        Check-ProjectInputs | Out-Null
        Invoke-WslStep
        Say ''
        Say 'Dry run: nothing was changed or started.'
        Done 0
    }

    Ensure-Docker
    $name = Get-ContainerName $script:ProjectDir
    if (Test-Running $name) {
        Fail "this project's analysis is already running. Watch it with:  run_atacseq logs   (or stop it: run_atacseq stop)"
    }
    # Work volume: an earlier run's work folder (metadata\work) is not reused by the volume.
    $onVolume = Get-ConfigBool 'work_volume' $true
    $volName = Get-WorkVolumeName $script:ProjectDir
    if ($onVolume) {
        $oldWork = Join-Path $script:ProjectDir 'metadata\work'
        $oldHas = (Test-Path -LiteralPath $oldWork) -and (@(Get-ChildItem -LiteralPath $oldWork -Force -ErrorAction SilentlyContinue).Count -gt 0)
        if ($oldHas -and -not (Test-WorkVolume $volName)) {
            Say ''
            Warn ("this project has an earlier work folder, metadata\work. With work_volume: true the run starts from " +
                  "scratch, because the resume state is in that folder. To resume it, set work_volume: false in " +
                  "$($script:ConfigName). To continue with the Docker work volume, answer yes (the old folder is " +
                  "left alone; run_atacseq clean deletes it).")
            if (-not (Confirm 'Start a fresh run with the Docker work volume?')) {
                Fail "not started. Set work_volume: false in $($script:ConfigName) to resume the earlier run."
            }
        } elseif ($oldHas) {
            Warn 'metadata\work from an earlier run is not used any more (the Docker work volume is). run_atacseq clean deletes it.'
        }
    }
    Repair-LineEndings
    $letters = @(Check-ProjectInputs)
    if ((Get-DriveMediaType $script:ProjectDir) -eq 'HDD') {
        if ($onVolume) {
            Warn ("this project folder is on a hard disk drive (HDD). The work files are on Docker's disk, so only " +
                  "the results are written here.")
        } else {
            Warn ("this project folder is on a hard disk drive (HDD). The analysis works, but it will be several " +
                  "times slower than on an SSD, and Windows may respond slowly meanwhile. An internal SSD is better if it has room.")
        }
    }

    # ATAC_HDD (conf/docker_tiers.config): 1 when a disk the run uses is a hard disk; the disk-heavy steps then run 2 at a time.
    $hdd = Get-RunHddValue $onVolume
    if ($hdd -eq '1') { Say '  A hard disk (HDD) holds the project, the inputs or the work files: the disk-heavy steps run 2 at a time (ATAC_HDD=1).' }

    Invoke-WslStep
    try {
        $budget = Get-RunBudget
        Ensure-Image

        $ts = Get-Date -Format 'yyyyMMdd_HHmmss'
        $imageId = "$(& docker image inspect $Image --format '{{.Id}}' 2>$null)".Trim()
        $configRel = $script:ConfigName
        # Mark the run as started now: until the container rewrites it, `status`
        # (and an agent polling it) would otherwise still see the previous run's outcome.
        $logDir = Join-Path $script:ProjectDir 'metadata\nf_logs'
        New-Item -ItemType Directory -Force $logDir | Out-Null
        Write-TextLF (Join-Path $logDir 'last_run.status') "running`n"

        # The work volume (idempotent: an existing one is kept as it is).
        if ($onVolume) {
            $null = & docker volume create --label "atacseq.project=$($script:ProjectDir)" $volName 2>&1
            if ($LASTEXITCODE -ne 0) {
                Write-TextLF (Join-Path $logDir 'last_run.status') "failed to start`n"
                Fail "Docker could not create the work volume $volName. Check that Docker is running, then run this again."
            }
        }

        $dockerArgs = @(
            'run', '-d', '--rm', '--init', '--name', $name, '--platform', 'linux/amd64',
            '--label', "atacseq.project=$($script:ProjectDir)",
            '--memory', "$($budget.memMB)m", '--memory-swap', "$($budget.memMB)m",
            '--cpuset-cpus', "0-$($budget.cpus - 1)",
            '-e', "ATAC_BUDGET_CPUS=$($budget.cpus)",
            '-e', "ATAC_BUDGET_MEMORY_MB=$($budget.memMB)",
            '-e', "ATAC_BUDGET_REASON=$($budget.reason)",
            '--mount', "type=bind,source=$($script:ProjectDir),target=/project",
            '--mount', "type=bind,source=$PipelineDir,target=/pipeline,readonly"
        )
        $dockerArgs += @(Get-MountArgs $letters)
        $dockerArgs += @(
            '-e', "ATAC_MOUNTS_HOST=$(Get-MountsHostText $letters)",
            '-e', "ATAC_IMAGE=$Image",
            '-e', "ATAC_IMAGE_ID=$imageId",
            '-e', "ATAC_PROJECT_HOST_DIR=$($script:ProjectDir)"
        )
        # Nextflow's work folder and resume cache: the project's Docker volume (/nfwork), or metadata\ (off).
        if ($onVolume) {
            $dockerArgs += @('--mount', "type=volume,source=$volName,target=/nfwork",
                '-e', 'ATAC_WORK_VOLUME=1', '-e', "ATAC_WORK_VOLUME_NAME=$volName")
        } else {
            $dockerArgs += @('-e', 'ATAC_WORK_VOLUME=0')
        }
        $dockerArgs += @('-e', "ATAC_HDD=$hdd")
        # Optional test settings, passed on when they are set on this computer.
        if ($env:ATAC_NF_PROFILE) { $dockerArgs += @('-e', "ATAC_NF_PROFILE=$($env:ATAC_NF_PROFILE)") }
        if ($env:ATAC_NF_EXTRA_ARGS) { $dockerArgs += @('-e', "ATAC_NF_EXTRA_ARGS=$($env:ATAC_NF_EXTRA_ARGS)") }
        if ($env:ATAC_WORK_DIR) { $dockerArgs += @('-e', "ATAC_WORK_DIR=$($env:ATAC_WORK_DIR)") }
        $dockerArgs += @($Image, 'bash', '/pipeline/docker/run_in_container.sh', $configRel, $ts)
        $null = & docker @dockerArgs
        if ($LASTEXITCODE -ne 0) {
            Write-TextLF (Join-Path $logDir 'last_run.status') "failed to start`n"
            Fail 'Docker could not start the pipeline container (see the message above).'
        }
        Start-Guard -container $name

        $log = "metadata\nf_logs\nextflow_$ts.log"
        Say ''
        Info 'The analysis is running inside Docker.'
        Say  "  Log file : $($script:ProjectDir)\$log"
        Say  "  Uses at most $($budget.cpus) CPU cores and $(Format-GB $budget.memMB) memory ($($budget.reason))."
        if ($onVolume) { Say  "  Work folder: Docker volume $volName (Docker's disk); results are written to $($script:ProjectDir)\results\" }
        else { Say  '  Work folder: metadata\work (project folder)' }
        Say  '  Progress : run_atacseq status      Stop: run_atacseq stop'
        Say  '  Closing this window does NOT stop the analysis. Keep the computer awake until it finishes.'
        if ($Detach) { Done 0 }
        Say  '  (Press Ctrl+C to stop watching; the analysis keeps running.)'
        Say ''
        Follow-Run $name
    } finally {
        # Also when the window is closed or Ctrl+C is pressed while following the log.
        if ($script:WslRestorePending) { $script:WslRestorePending = $false; Restore-WslIfIdle }
    }
}

function Cmd-Status {
    Get-ProjectDirFromArgs
    $name = Get-ContainerName $script:ProjectDir
    $running = $false
    if (Test-DockerUp) { $running = Test-Running $name }
    $statusFile = Join-Path $script:ProjectDir 'metadata\nf_logs\last_run.status'
    $logs = @(Get-ChildItem (Join-Path $script:ProjectDir 'metadata\nf_logs') -Filter 'nextflow_*.log' -ErrorAction SilentlyContinue | Sort-Object Name)
    if ($logs.Count -eq 0) { Say "No run yet in $($script:ProjectDir)."; Show-WorkFolderStatus; Done 0 }
    $last = $logs[-1]
    if ($running) { Info 'Status: RUNNING' }
    elseif (Test-Path $statusFile) {
        $st = ([System.IO.File]::ReadAllText($statusFile)).Trim()
        if ($st -eq 'running') { $st = 'interrupted (Docker or the computer stopped mid-run; run again to resume)' }
        Info "Status: last run $st"
    }
    Say "Log: $($last.FullName)"
    Show-WorkFolderStatus

    # The trace is named by the launch timestamp (nextflow.config): the newest one is this project's latest run.
    $traces = @(Get-ChildItem (Join-Path $script:ProjectDir 'metadata\pipeline_info') -Filter 'execution_trace*.txt' -ErrorAction SilentlyContinue | Sort-Object Name)
    $trace = if ($traces.Count -gt 0) { $traces[-1].FullName } else { '' }
    if ($trace -and (Test-Path $trace)) {
        $rows = @(Import-Csv $trace -Delimiter "`t")
        if ($rows.Count -gt 0) {
            Say ''
            Say 'Steps (this run):'
            $rows | Group-Object status | ForEach-Object { Say ("  {0,-10} {1}" -f $_.Name, $_.Count) }
            $failed = @($rows | Where-Object { $_.status -eq 'FAILED' })
            foreach ($f in $failed) { Say "  FAILED: $($f.name) (exit $($f.exit))" }
        }
    }
    $ft = Join-Path $script:ProjectDir ('metadata\failed_tasks\' + ($last.BaseName -replace '^nextflow_', ''))
    if (Test-Path $ft) { Say ''; Say "Files of the failed step(s): $ft" }
    Say ''
    Say "Last lines of the log:"
    Get-Content $last.FullName -Tail 25 | ForEach-Object { Say "  $_" }
    Done 0
}

function Cmd-Logs {
    Get-ProjectDirFromArgs
    $name = Get-ContainerName $script:ProjectDir
    Ensure-Docker
    if (Test-Running $name) { Say '(Ctrl+C stops watching; the analysis keeps running.)'; Follow-Run $name }
    $logs = @(Get-ChildItem (Join-Path $script:ProjectDir 'metadata\nf_logs') -Filter 'nextflow_*.log' -ErrorAction SilentlyContinue | Sort-Object Name)
    if ($logs.Count -eq 0) { Say 'No run yet.'; Done 0 }
    Say "Not running. Last log: $($logs[-1].FullName)"
    Get-Content $logs[-1].FullName -Tail 60
    Done 0
}

function Cmd-Stop {
    Get-ProjectDirFromArgs
    $name = Get-ContainerName $script:ProjectDir
    Ensure-Docker
    if (-not (Test-Running $name)) { Say 'Nothing is running for this project.'; Done 0 }
    Info 'Stopping the analysis (finished steps are kept; running it again resumes)...'
    $null = & docker stop -t 60 $name 2>&1
    $statusFile = Join-Path $script:ProjectDir 'metadata\nf_logs\last_run.status'
    Write-TextLF $statusFile "stopped by the user`n"
    Say 'Stopped.'
    Done 0
}

function Cmd-Clean {
    Get-ProjectDirFromArgs
    $name = Get-ContainerName $script:ProjectDir
    # metadata\ holds the work folder and the resume cache; the project's Docker work volume holds them too (work_volume: true).
    $volName = Get-WorkVolumeName $script:ProjectDir
    $targets = @(
        (Join-Path $script:ProjectDir 'metadata\work'),
        (Join-Path $script:ProjectDir 'metadata\.nextflow')
    )
    $existing = @($targets | Where-Object { Test-Path -LiteralPath $_ })
    if ($existing.Count -eq 0 -and -not (Test-DockerUp)) {
        Say 'Nothing to clean in metadata\ (no work folder or resume cache).'
        Say "The Docker work volume ($volName) is skipped: Docker is not running. Start Docker and run clean again to delete it."
        Done 0
    }
    Ensure-Docker
    $volExists = Test-WorkVolume $volName
    if ($existing.Count -eq 0 -and -not $volExists) {
        Say 'Nothing to clean (no metadata\work, metadata\.nextflow folder or Docker work volume).'
        Done 0
    }
    if (Test-Running $name) { Fail 'the analysis is running - stop it first (run_atacseq stop).' }
    $bytes = [double]0
    foreach ($t in $existing) {
        $sum = (Get-ChildItem -LiteralPath $t -Recurse -Force -File -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
        if ($sum) { $bytes += [double]$sum }
    }
    if ($volExists) {
        $volBytes = Get-WorkVolumeBytes $volName
        if ($null -ne $volBytes) { $bytes += [double]$volBytes }
    }
    $parts = @($existing | ForEach-Object { "$_" })
    if ($volExists) { $parts += "the Docker work volume $volName" }
    $shown = $parts -join ', '
    Say "This deletes ${shown}: Nextflow's working copies of every step and its resume cache, about $(Format-Bytes $bytes)."
    Say 'results\, original_data\ and your config and CSV files are NOT touched, but the next run will start from scratch.'
    if (-not (Confirm 'Delete it?')) { Say 'Nothing deleted.'; Done 0 }
    if ((Quiet { docker image inspect $Image }) -eq 0) {
        # Deleted from inside Linux: the work folder holds Linux symlinks that Windows tools handle poorly.
        $null = & docker run --rm --platform linux/amd64 --mount "type=bind,source=$($script:ProjectDir),target=/project" $Image rm -rf /project/metadata/work /project/metadata/.nextflow 2>&1
    }
    foreach ($t in $existing) {
        if (Test-Path -LiteralPath $t) { Remove-Item -LiteralPath $t -Recurse -Force -ErrorAction SilentlyContinue }
    }
    if ($volExists) {
        if ((Quiet { docker volume rm $volName }) -eq 0) { Say "Deleted the Docker work volume $volName." }
        else { Warn "could not delete the Docker work volume $volName (is a container still using it?)" }
    }
    Say "Deleted. About $(Format-Bytes $bytes) of disk space is free again."
    Done 0
}

# run_atacseq script utils\X.R [ARGS]: one of the project's own scripts, run in the
# pipeline's software environment. Foreground only; the project is the current folder.
function Cmd-Script {
    $script:ProjectDir = (Get-Location).ProviderPath.TrimEnd('\')
    $script:ConfigName = 'config.yaml'
    $script:ConfigFile = Join-Path $script:ProjectDir $script:ConfigName
    if ($positional.Count -eq 0 -or "$($positional[0])".Trim() -eq '') {
        Fail 'usage: run_atacseq script utils\<name>.<R|py|sh> [args...]   (run from the project folder)'
    }
    $given = [string]$positional[0]
    $candidate = $given
    if (-not [System.IO.Path]::IsPathRooted($given)) { $candidate = Join-Path $script:ProjectDir ($given -replace '/', '\') }
    if (-not (Test-Path -LiteralPath $candidate)) { Fail "$given not found (looked for $candidate)." }
    if (-not (Test-Path -LiteralPath $candidate -PathType Leaf)) { Fail "$given is not a file." }
    $full = (Resolve-Path -LiteralPath $candidate).ProviderPath
    $root = $script:ProjectDir + '\'
    if (-not $full.StartsWith($root, [StringComparison]::OrdinalIgnoreCase)) {
        Fail "The script must be inside the project folder ($($script:ProjectDir))."
    }
    if (([System.IO.Path]::GetExtension($full)).ToLowerInvariant() -notin @('.r', '.py', '.sh')) {
        Fail 'unsupported extension (use .R, .py or .sh)'
    }
    $rel = Get-ProjectRelative $full

    Ensure-Docker
    Repair-LineEndings
    $letters = @(Check-InputDrives)
    $letters = @(Add-ArgDrives $letters @($positional | Select-Object -Skip 1))
    $budget = Get-RunBudget
    Ensure-Image
    if (Test-Running (Get-ContainerName $script:ProjectDir)) {
        Warn "the analysis of this project is running and shares this computer's memory with the script (the script may run slower)."
    }

    $ts = Get-Date -Format 'yyyyMMdd_HHmmss'
    $base = [System.IO.Path]::GetFileNameWithoutExtension($full)
    $logName = "$($base)_$($ts).log"
    $logRel = "metadata/script_logs/$logName"
    $logDir = Join-Path $script:ProjectDir 'metadata\script_logs'
    New-Item -ItemType Directory -Force $logDir | Out-Null
    $logWin = Join-Path $logDir $logName
    $imageId = "$(& docker image inspect $Image --format '{{.Id}}' 2>$null)".Trim()

    Info "Running $rel in the pipeline's software environment (at most $($budget.cpus) CPU cores and $(Format-GB $budget.memMB) memory). Log: $logWin"
    $dockerArgs = @(
        'run', '--rm', '--init', '--platform', 'linux/amd64',
        '--label', "atacseq.project=$($script:ProjectDir)",
        '--memory', "$($budget.memMB)m", '--memory-swap', "$($budget.memMB)m",
        '--cpuset-cpus', "0-$($budget.cpus - 1)",
        '-e', "ATAC_BUDGET_CPUS=$($budget.cpus)",
        '-e', "ATAC_BUDGET_MEMORY_MB=$($budget.memMB)",
        '-e', "ATAC_BUDGET_REASON=$($budget.reason)",
        '-e', "ATAC_IMAGE=$Image",
        '-e', "ATAC_IMAGE_ID=$imageId",
        '-e', "ATAC_PROJECT_HOST_DIR=$($script:ProjectDir)",
        '-e', "ATAC_MOUNTS_HOST=$(Get-MountsHostText $letters)",
        '--mount', "type=bind,source=$($script:ProjectDir),target=/project",
        '--mount', "type=bind,source=$PipelineDir,target=/pipeline,readonly"
    )
    $dockerArgs += @(Get-MountArgs $letters)
    $dockerArgs += @('-w', '/project', $Image, 'bash', '/pipeline/docker/run_script_in_container.sh', $rel, $logRel)
    $dockerArgs += @($positional | Select-Object -Skip 1)
    # Windows PowerShell 5.1 strips double quotes from native-command arguments, so script arguments cannot contain them.
    & docker @dockerArgs
    $rc = $LASTEXITCODE

    Say ''
    $msg = "Script finished (exit code $rc). Log: $logWin"
    if ($rc -eq 0) { Info $msg } else { Write-Host $msg -ForegroundColor Red }
    Done $rc
}

# Hidden (maintainers): the Docker memory plan for hypothetical inputs, without Docker.
#   wsl-plan [-TotalRamGB N] [-CurrentGB N|none|nofile]
function Cmd-WslPlan {
    $ram = $null; $curArg = 'none'
    for ($i = 0; $i -lt $positional.Count; $i++) {
        if ($positional[$i] -match '^-+TotalRamGB$' -and ($i + 1) -lt $positional.Count) {
            $ram = [double]::Parse([string]$positional[$i + 1], [System.Globalization.CultureInfo]::InvariantCulture); $i++
        } elseif ($positional[$i] -match '^-+CurrentGB$' -and ($i + 1) -lt $positional.Count) {
            $curArg = [string]$positional[$i + 1]; $i++
        } else {
            Fail "wsl-plan: unknown argument $($positional[$i])"
        }
    }
    if ($null -eq $ram) { $ram = (Get-HostResources).memMB / 1024.0; Say "(no -TotalRamGB given: this computer's memory is used)" }
    $tail = "swap=2GB`r`n`r`n[experimental]`r`nautoMemoryReclaim=gradual`r`n"
    if ($curArg -eq 'nofile') { $text = '' }
    elseif ($curArg -eq 'none') { $text = "[wsl2]`r`n$tail" }
    elseif ($curArg -match '^\d+(\.\d+)?$') { $text = "[wsl2]`r`nmemory=${curArg}GB`r`n$tail" }
    else { Fail "wsl-plan: -CurrentGB takes a number of GB, 'none' (no memory= line) or 'nofile' (no .wslconfig)." }
    Say ("Plan for a computer with {0:0.#} GB of memory. Nothing is read from or written to Docker or .wslconfig." -f $ram)
    $plan = Get-WslPlan $text $ram 0 $false $false
    Show-WslPlanReport $plan $ram $text $false
    Done 0
}

# Hidden (maintainers): the drive letters a run would mount; refuses network paths.
function Cmd-ScanPaths {
    Get-ProjectDirFromArgs
    if (-not (Test-Path -LiteralPath $script:ConfigFile)) { Fail "$($script:ConfigName) not found in $($script:ProjectDir)." }
    Say "Scanning $($script:ConfigName) and the CSV files in $($script:ProjectDir)"
    $letters = @(Check-InputDrives)
    if ($letters.Count -eq 0) { Say '  No Windows drive is named: only the project folder is visible inside the analysis.' }
    foreach ($l in $letters) { Say "  would mount ${l}:\ read-only as /mnt/$($l.ToLowerInvariant())" }
    Check-OutputDirs
    Check-ConfigFiles
    Done 0
}

# Hidden (maintainers): runs the pure .wslconfig functions on sample files.
function Cmd-SelfTest {
    $crlf = "`r`n"
    $cases = @(
        @{ Name = 'no .wslconfig yet';                       Text = '' },
        @{ Name = '[wsl2] without memory (LF)';              Text = "[wsl2]`nswap=2GB`n`n[experimental]`nautoMemoryReclaim=gradual`n" },
        @{ Name = 'memory with comments and spacing (CRLF)'; Text = "# my Docker settings${crlf}[wsl2]${crlf}; keep the swap line${crlf}memory = 6GB  # six gigabytes${crlf}swap=2GB${crlf}" },
        @{ Name = 'memory= only in another section';         Text = "[experimental]`nmemory=1GB`n`n[wsl2]`nswap=2GB`n" },
        @{ Name = 'typical laptop file (LF)';                Text = "[wsl2]`nmemory=6GB`nswap=2GB`n`n[experimental]`nautoMemoryReclaim=gradual`n" },
        @{ Name = 'memory line only (CRLF)';                 Text = "[wsl2]${crlf}memory=6GB${crlf}swap=2GB${crlf}" },
        @{ Name = 'no trailing newline';                     Text = "[wsl2]`nmemory=6GB`nswap=2GB" },
        @{ Name = 'no [wsl2] section, no trailing newline';  Text = "[experimental]`nautoMemoryReclaim=gradual" },
        @{ Name = 'no memory= line (CRLF), no trailing newline'; Text = "[wsl2]${crlf}swap=2GB${crlf}${crlf}[experimental]${crlf}autoMemoryReclaim=gradual" }
    )
    $failed = 0
    foreach ($c in $cases) {
        $orig = [string]$c.Text
        $new = Add-WslMemoryText $orig 16
        $problems = @()
        if ((Get-WslMemoryValue $new) -ne '16GB') { $problems += "memory value is '$(Get-WslMemoryValue $new)', not 16GB" }
        if ((Find-MarkerIndex @(Split-TextLines $new)) -lt 0) { $problems += 'no marker line' }
        if ((Remove-WslMemoryText $new).Text -cne $orig) { $problems += 'restore is not byte-identical to the original' }
        if ($orig.Contains($crlf) -and ((($new -replace $crlf, '')).Contains("`n"))) { $problems += 'a bare LF was written into a CRLF file' }
        if ($problems.Count -eq 0) { Say "PASS  $($c.Name)" }
        else { $failed++; Write-Host "FAIL  $($c.Name): $($problems -join '; ')" -ForegroundColor Red }
    }
    # Restore must leave a file without a marker alone.
    $plain = "[wsl2]`nmemory=6GB`n"
    if ((Remove-WslMemoryText $plain).Text -ceq $plain) { Say 'PASS  restore on a file without a marker changes nothing' }
    else { $failed++; Write-Host 'FAIL  restore on a file without a marker changed it' -ForegroundColor Red }
    # Docker work volume helpers (pure: no Docker needed).
    $vol = Get-WorkVolumeName 'D:\Projects\ATAC_X'
    $volChecks = @(
        @{ Name = 'work volume name is deterministic';                 Ok = ((Get-WorkVolumeName 'D:\Projects\ATAC_X') -eq $vol) },
        @{ Name = 'same folder in another case gives the same name';   Ok = ((Get-WorkVolumeName 'd:\projects\atac_x') -eq $vol) },
        @{ Name = 'work volume name is a valid Docker name';           Ok = ($vol -match '^[a-zA-Z0-9][a-zA-Z0-9_.-]*$') },
        @{ Name = 'work volume name contains the folder name, ends _work'; Ok = (($vol -like '*atac_x*') -and $vol.EndsWith('_work')) },
        @{ Name = 'publish symlink with work volume is refused';       Ok = ((Get-PublishModeProblem 'symlink' $true) -ne '') },
        @{ Name = 'publish rellink with work volume is refused';       Ok = ((Get-PublishModeProblem 'rellink' $true) -ne '') },
        @{ Name = 'publish link with work volume is refused';          Ok = ((Get-PublishModeProblem 'link' $true) -ne '') },
        @{ Name = 'publish auto with work volume is allowed';          Ok = ((Get-PublishModeProblem 'auto' $true) -eq '') },
        @{ Name = 'publish copy with work volume is allowed';          Ok = ((Get-PublishModeProblem 'copy' $true) -eq '') },
        @{ Name = 'publish symlink without work volume is allowed';    Ok = ((Get-PublishModeProblem 'symlink' $false) -eq '') },
        @{ Name = 'df: Available column of a 2-line df -Pk';           Ok = ((Get-DfAvailBytes ("Filesystem 1024-blocks Used Available Capacity Mounted on`noverlay 1055762868 20000000 982000000 3% /")) -eq ([double]982000000 * 1024)) },
        @{ Name = 'df: garbage gives nothing';                         Ok = ($null -eq (Get-DfAvailBytes 'hello there')) },
        @{ Name = 'Docker data folder from DataFolder';                Ok = ((Get-DockerDataFolder '{"DataFolder":"D:\\DockerData"}') -eq 'D:\DockerData') },
        @{ Name = 'Docker data folder from dataFolder';                Ok = ((Get-DockerDataFolder '{"dataFolder":"E:\\x"}') -eq 'E:\x') },
        @{ Name = 'Docker data folder: invalid JSON gives nothing';    Ok = ($null -eq (Get-DockerDataFolder '{ not json')) },
        @{ Name = 'disk fit: room to spare gives no message';          Ok = ((Get-DiskFit (100*1GB) (50*1GB)) -eq '') },
        @{ Name = 'disk fit: 10% to spare is a Warn';                  Ok = ((Get-DiskFit (55*1GB) (50*1GB)) -eq 'warn') },
        @{ Name = 'disk fit: exactly the need is a Warn';              Ok = ((Get-DiskFit (50*1GB) (50*1GB)) -eq 'warn') },
        @{ Name = 'disk fit: below the need is a Fail';                Ok = ((Get-DiskFit (40*1GB) (50*1GB)) -eq 'fail') },
        @{ Name = 'disk fit: nothing needed gives no message';         Ok = ((Get-DiskFit 0 0) -eq '') },
        @{ Name = 'ATAC_HDD: one HDD gives 1';                         Ok = ((Get-AtacHddValue @('SSD', 'HDD')) -eq '1') },
        @{ Name = 'ATAC_HDD: SSD and unknown give 0';                  Ok = ((Get-AtacHddValue @('SSD', '')) -eq '0') },
        @{ Name = 'ATAC_HDD: no drives known gives 0';                 Ok = ((Get-AtacHddValue @()) -eq '0') }
    )
    foreach ($c in $volChecks) {
        if ($c.Ok) { Say "PASS  $($c.Name)" }
        else { $failed++; Write-Host "FAIL  $($c.Name)" -ForegroundColor Red }
    }
    # Sample plan: the typical laptop file at 16 GB of RAM (tier 8 GB).
    Say ''
    Say 'Example (typical laptop file, 16 GB of RAM):'
    $example = "[wsl2]`nmemory=6GB`nswap=2GB`n`n[experimental]`nautoMemoryReclaim=gradual`n"
    Show-WslDiff $example (Add-WslMemoryText $example (Get-WslTierGB 16))
    Say ''
    if ($failed -eq 0) { Info 'All self-test cases passed.'; Done 0 }
    Write-Host "$failed self-test case(s) FAILED." -ForegroundColor Red
    Done 1
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
if ($Command -match '^(setup|run|script|clean|status)$') { $Image = Get-ImageTag }
# A run that was closed before Docker's memory limit was put back: offer it now.
if ($Command -notmatch '^(help|-h|--help|/\?|selftest|wsl-plan|scan-paths)$') { Restore-WslIfIdle -Ask }

switch -Regex ($Command) {
    '^setup$'       { Cmd-Setup }
    '^init$'        { Cmd-Init }
    '^run$'         { Cmd-Run }
    '^status$'      { Cmd-Status }
    '^logs$'        { Cmd-Logs }
    '^stop$'        { Cmd-Stop }
    '^clean$'       { Cmd-Clean }
    '^script$'      { Cmd-Script }
    '^wsl-plan$'    { Cmd-WslPlan }
    '^scan-paths$'  { Cmd-ScanPaths }
    '^selftest$'    { Cmd-SelfTest }
    default         { Show-Help; Done 0 }
}
