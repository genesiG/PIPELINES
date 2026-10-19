# run_crisprscreen.ps1 - CRISPR_SCREEN pipeline launcher, Docker version for Windows.
#
# Labmates (and AI agents) call it through run_crisprscreen.bat, which only
# exists so PowerShell's execution policy never gets in the way:
#
#   run_crisprscreen setup            one-time: check/install Docker Desktop, build the
#                                     pipeline's Docker image, put this command on PATH
#   run_crisprscreen init [FOLDER]    set up a new project folder (default: the current folder)
#   run_crisprscreen [run] [CONFIG]   run the pipeline for the project in the current folder
#                                     (CONFIG defaults to config.yaml)
#   run_crisprscreen status           is it running? how far along? did it fail, and where?
#   run_crisprscreen logs             follow the live log of the running analysis
#   run_crisprscreen stop             stop a running analysis (re-running resumes it)
#   run_crisprscreen clean            delete work/ now (normally removed after a successful run unless keep_work: true)
#   run_crisprscreen help
#
# Options: --detach (run: start and return, don't follow the log), --yes
# (answer yes to questions: for AI agents / scripts), --pause (wait for Enter
# before closing; used by the double-clickable RUN_ANALYSIS.bat), --rebuild
# (setup: rebuild the Docker image even if it is up to date), --docker-memory
# (setup: let Docker Desktop use the full memory budget, editing .wslconfig and
# restarting WSL; --yes alone never does this, it only says how).
#
# What `run` does (the Docker counterpart of the cluster's run_crisprscreen.sh):
#   1. Checks the raw FASTQ path(s) named in config.yaml. The FASTQs stay where
#      they are: they are never moved, copied or renamed, and config.yaml is not
#      rewritten for them. Outside the project folder they are bind-mounted
#      read-only under /inputs/... (UNC paths and network drives are refused).
#      Any other path in config.yaml must be inside the project folder (it is
#      rewritten as a relative path if needed).
#   2. Makes sure Docker Desktop is running and the image is built and up to
#      date with docker/ (rebuilt automatically after a pipeline update).
#   3. Starts one container that runs Nextflow and every step
#      (docker/run_in_container.sh), with the project folder mounted at
#      /project and this pipeline folder (read-only) at /pipeline, then
#      follows its log. Closing the window does NOT stop the analysis.

$ErrorActionPreference = 'Continue'
Set-StrictMode -Version 2

$PipelineDir  = $PSScriptRoot
$PipelineName = 'crisprscreen'
$Image        = 'rephilllab/crisprscreen:0.1.0'
if ($env:CRISPR_IMAGE_OVERRIDE) { $Image = $env:CRISPR_IMAGE_OVERRIDE }

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
$script:Pause = $false; $script:Yes = $false; $Detach = $false; $Rebuild = $false; $DockerMemory = $false
$positional = @()
foreach ($a in $args) {
    switch -Regex ([string]$a) {
        '^--?pause$'   { $script:Pause = $true; continue }
        '^--?docker-memory$' { $DockerMemory = $true; continue }
        '^--?yes$|^-y$' { $script:Yes = $true; continue }
        '^--?detach$'  { $Detach = $true; continue }
        '^--?rebuild$' { $Rebuild = $true; continue }
        default        { $positional += [string]$a }
    }
}
$Command = 'run'
if ($positional.Count -gt 0) {
    if ($positional[0] -match '^(setup|init|run|status|logs|stop|clean|help|-h|--help|/\?)$') {
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
              "(choose 'Use WSL 2'), restart the computer if asked, start Docker Desktop once, then run this again.")
    }
    if (-not (Confirm 'Install Docker Desktop now with winget? (Windows will ask for administrator permission)')) {
        Fail 'Docker Desktop is required. Install it from https://www.docker.com/products/docker-desktop/ and run this again.'
    }
    Info 'Installing Windows Subsystem for Linux (WSL 2) and Docker Desktop...'
    & wsl.exe --install --no-distribution
    & winget install -e --id Docker.DockerDesktop --accept-package-agreements --accept-source-agreements
    if ($LASTEXITCODE -ne 0) { Fail "winget could not install Docker Desktop (exit $LASTEXITCODE)." }
    Say ''
    Say 'Docker Desktop is installed. Next:'
    Say '  1. RESTART the computer (needed once, for WSL 2).'
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
        Fail "Docker is not installed. Run:  run_crisprscreen setup"
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

# SHA-256 of what the image is built from (docker/Dockerfile + docker/env/*.explicit.txt),
# line endings normalized; stored as an image label, so a changed docker/
# folder (e.g. after `git pull`) triggers a rebuild.
function Get-DockerContextHash {
    $files = @(Get-Item (Join-Path $PipelineDir 'docker\Dockerfile'))
    $files += @(Get-ChildItem (Join-Path $PipelineDir 'docker\env') -Filter '*.explicit.txt' | Sort-Object Name)
    $ms = New-Object System.IO.MemoryStream
    foreach ($f in $files) {
        $name = [System.Text.Encoding]::UTF8.GetBytes($f.Name + "`n")
        $ms.Write($name, 0, $name.Length)
        $bytes = [System.IO.File]::ReadAllBytes($f.FullName) | Where-Object { $_ -ne 13 }
        if ($bytes) { $arr = [byte[]]$bytes; $ms.Write($arr, 0, $arr.Length) }
    }
    $sha = [System.Security.Cryptography.SHA256]::Create()
    $hash = $sha.ComputeHash($ms.ToArray())
    return (($hash | ForEach-Object { $_.ToString('x2') }) -join '')
}

function Ensure-Image([switch]$Force) {
    $want = Get-DockerContextHash
    # Labels as JSON: Windows PowerShell 5.1 strips double quotes inside native
    # arguments, so a template like {{ index .Config.Labels "x" }} never reaches docker intact.
    $labels = (& docker image inspect $Image --format '{{json .Config.Labels}}' 2>$null)
    $exists = ($LASTEXITCODE -eq 0)
    $have = ''
    if ($exists) { try { $have = "$(($labels | ConvertFrom-Json).'crisprscreen.docker_hash')".Trim() } catch { $have = '' } }
    if (-not $Force -and $exists -and $have -eq $want) { return }
    if ($exists -and -not $Force) {
        Info "The pipeline's Docker image is out of date (docker/ changed) - rebuilding it."
    } else {
        Info "Building the pipeline's Docker image ($Image)."
        Say  'This happens once (and again only after a pipeline update changes docker/).'
        Say  'It downloads about 3 GB and takes roughly 10-30 minutes. Keep the computer awake and online.'
    }
    $ctx = Join-Path $PipelineDir 'docker'
    Start-Guard -watchPid $PID
    & docker build --platform linux/amd64 --progress=plain --label "crisprscreen.docker_hash=$want" -t $Image $ctx
    $rc = $LASTEXITCODE
    # The build leaves gigabytes of image layers cached in Docker's VM: hand that memory back now.
    $null = & wsl.exe -d docker-desktop -u root sh -c 'echo 1 > /proc/sys/vm/drop_caches' 2>&1
    if ($rc -ne 0) {
        Fail ("Building the Docker image failed (see the messages above). Check the internet connection " +
              "and free disk space, then run this again.")
    }
    Info "Image $Image is ready."
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
# config.yaml: flat keys only (as on the cluster), so simple line edits do.
# ---------------------------------------------------------------------------
function Get-YamlValue([string]$file, [string]$key) {
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

function Set-YamlValue([string]$file, [string]$key, [string]$value) {
    $text = [System.IO.File]::ReadAllText($file)
    $nl = "`n"; if ($text.Contains("`r`n")) { $nl = "`r`n" }
    $lines = $text -split "`r?`n"
    $quoted = $value
    if ($value -match '[:#\s]' -or $value -match '^[\[\{&*!|>%@`]') { $quoted = "'" + ($value -replace "'", "''") + "'" }
    for ($i = 0; $i -lt $lines.Count; $i++) {
        if ($lines[$i] -match "^$([regex]::Escape($key)):") { $lines[$i] = "${key}: $quoted"; break }
    }
    Write-TextLF $file ($lines -join $nl)
}

function Test-Unset([string]$v) { return (-not $v) -or ($v -eq 'null') -or ($v -eq '~') }

function Test-AbsolutePath([string]$v) {
    return ($v -match '^[A-Za-z]:[\\/]') -or ($v -match '^\\\\') -or ($v -match '^/') -or ($v -match '^~')
}

# A path from config.yaml, as a full Windows path (relative = relative to the project).
function Resolve-ConfigPath([string]$v) {
    if ($v -match '^~[\\/]?(.*)$') { $v = Join-Path $HOME $Matches[1] }
    elseif ($v -match '^/([A-Za-z])/(.*)$') { $v = "$($Matches[1]):\$($Matches[2])" }   # Git Bash style /d/data/...
    elseif ($v -match '^/mnt/([A-Za-z])/(.*)$') { $v = "$($Matches[1]):\$($Matches[2])" } # WSL style
    if ($v -match '^[A-Za-z]:[\\/]' -or $v -match '^\\\\') { return [System.IO.Path]::GetFullPath(($v -replace '/', '\')) }
    return [System.IO.Path]::GetFullPath((Join-Path $script:ProjectDir ($v -replace '/', '\')))
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

# A raw FASTQ path named in config.yaml: fastq_r1 / fastq_r2 are files, fastq_dir
# is a folder. The data is never moved, copied or renamed, and config.yaml is not
# rewritten for these keys. Returns $null if the key is not set.
function Resolve-RawFastq([string]$key) {
    $v = Get-YamlValue $script:ConfigFile $key
    if (Test-Unset $v) { return $null }
    $full = Resolve-ConfigPath $v
    if ($key -eq 'fastq_dir') {
        if (-not (Test-Path -LiteralPath $full -PathType Container)) {
            Fail "$key in $($script:ConfigName) points to a folder that doesn't exist: $v (looked for $full)"
        }
        $files = @(Get-ChildItem -LiteralPath $full -File | Where-Object { $_.Name -match '\.fastq(\.gz)?$' })
        if ($files.Count -eq 0) { Fail "fastq_dir $v contains no *.fastq.gz or *.fastq files (looked in $full)." }
    } elseif (-not (Test-Path -LiteralPath $full -PathType Leaf)) {
        Fail "$key in $($script:ConfigName) points to a file that doesn't exist: $v (looked for $full)"
    }
    $fq = @{ key = $key; full = $full; inside = $false; hostDir = ''; target = ''; containerPath = '' }
    $rel = Get-ProjectRelative $full
    if ($rel) {
        $fq.inside = $true
        if ($rel -eq '.') { $fq.containerPath = '/project' } else { $fq.containerPath = '/project/' + $rel }
        return $fq
    }

    # Outside the project folder the container cannot see it unless it is mounted.
    $network = $false
    try { $network = ([System.IO.DriveInfo]::new([System.IO.Path]::GetPathRoot($full)).DriveType -eq [System.IO.DriveType]::Network) } catch { }
    if ($full.StartsWith('\\') -or $network) {
        Fail "$key ($full) is on a network share. Docker Desktop cannot mount network shares: copy the file (or folder) to a local drive or a USB drive and point config.yaml to it."
    }
    if ($full.Contains(',')) {
        Fail "$key ($full) has a comma in its path, which Docker cannot mount. Rename or move it (remove the comma) and point config.yaml to the new path."
    }
    if ($key -eq 'fastq_dir') {
        $fq.hostDir = $full; $fq.target = '/inputs/fastq_dir'; $fq.containerPath = '/inputs/fastq_dir'
    } else {
        $fq.hostDir = [System.IO.Path]::GetDirectoryName($full)
        $fq.target = "/inputs/$key"
        $fq.containerPath = "/inputs/$key/" + [System.IO.Path]::GetFileName($full)
    }
    return $fq
}

# Docker arguments that make one FASTQ path visible inside the container: outside
# the project a read-only bind mount + two env vars; inside it only the env vars (/project/...).
function Get-FastqInputArgs($fq) {
    $out = @()
    $envName = 'CRISPR_' + $fq.key.ToUpperInvariant()
    if ($fq.inside) {
        $out += @('-e', "$envName=$($fq.containerPath)", '-e', "${envName}_HOST=$($fq.full)")
        Say "  $($fq.key): $($fq.full) (in the project folder)"
        return $out
    }
    $out += @('--mount', "type=bind,source=$($fq.hostDir),target=$($fq.target),readonly")
    $out += @('-e', "$envName=$($fq.containerPath)", '-e', "${envName}_HOST=$($fq.full)")
    Say "  $($fq.key): $($fq.full) (stays where it is; mounted read-only)"
    return $out
}

# Every other path in config.yaml must be inside the project folder (the
# container sees nothing else); written relative, with forward slashes.
function Check-ConfigPaths {
    $pathKeys = @('samplesheet', 'library', 'contrasts', 'outdir', 'original_data_dir', 'intermediate_dir', 'metadata_dir')
    foreach ($line in [System.IO.File]::ReadAllLines($script:ConfigFile)) {
        if ($line -match '^([A-Za-z0-9_]+):\s*(.+?)\s*(#.*)?$') {
            $k = $Matches[1]; $raw = $Matches[2].Trim("'`"")
            if ($k -in @('fastq_r1', 'fastq_r2', 'fastq_dir')) { continue }
            if (($k -in $pathKeys -or $raw -match '\.(r|R|py|sh)$') -and -not (Test-Unset $raw)) { $pathKeys += $k }
        }
    }
    foreach ($k in ($pathKeys | Select-Object -Unique)) {
        $v = Get-YamlValue $script:ConfigFile $k
        if (Test-Unset $v) { continue }
        if (Test-AbsolutePath $v) {
            $rel = Get-ProjectRelative (Resolve-ConfigPath $v)
            if (-not $rel) {
                Fail ("${k}: $v is outside the project folder, which the pipeline cannot see. Copy it into " +
                      "$($script:ProjectDir) and set ${k} to its name there (e.g. ${k}: $([System.IO.Path]::GetFileName($v)))")
            }
            Say "  ${k}: $v -> $rel (relative to the project folder)"
            Set-YamlValue $script:ConfigFile $k $rel
        } elseif ($v.Contains('\')) {
            Set-YamlValue $script:ConfigFile $k ($v -replace '\\', '/')
        }
    }
    foreach ($k in @('samplesheet')) {
        $v = Get-YamlValue $script:ConfigFile $k
        if (-not (Test-Unset $v) -and -not (Test-Path -LiteralPath (Resolve-ConfigPath $v))) {
            Fail "${k}: $v not found in the project folder. Create it (see examples/samplesheet_example.csv)."
        }
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

# ---------------------------------------------------------------------------
# Resources: one budget per run, enforced as the container's hard limits.
#
#   memory = 75% of the computer's RAM, but always leaving 3 GB to Windows
#            (a plain percentage starves small laptops: 80% of 8 GB leaves 1.6 GB)
#   CPUs   = all logical cores but 1 (but 2 from 8 cores up), so the desktop
#            stays responsive; heavy CPU use only slows a computer, it can't harm it
#   ... and never more than Docker Desktop itself may use (its WSL 2 VM).
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

function Format-GB([int]$mb) { return ('{0:0.#} GB' -f ($mb / 1024.0)) }

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
            Fail "max_memory_gb: $wantMem in $($script:ConfigName) is below the $($minMem / 1024) GB this pipeline needs. Raise it or remove the line."
        } elseif ($p.memMB -ge $minMem) {
            Fail ("Docker Desktop may only use $(Format-GB $d.memMB), and this pipeline needs about $($minMem / 1024) GB plus Docker's own share. " +
                  "Run  run_crisprscreen setup  - it offers to let Docker use more of this computer's memory.")
        } else {
            Fail ("This computer has $(Format-GB $h.memMB) of memory. Keeping 3 GB for Windows leaves less than the $($minMem / 1024) GB " +
                  "this pipeline needs. Run it on a computer with at least $([math]::Ceiling(($minMem + 3072) / 1024)) GB, " +
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

# .wslconfig [wsl2] memory=: how much memory Docker Desktop's WSL 2 VM may use.
function Set-WslMemory([int]$gb) {
    $f = Join-Path $env:USERPROFILE '.wslconfig'
    $lines = @()
    if (Test-Path $f) { $lines = @([System.IO.File]::ReadAllLines($f)) }
    $out = New-Object System.Collections.Generic.List[string]
    $inWsl2 = $false; $sawWsl2 = $false; $done = $false
    foreach ($l in $lines) {
        if ($l -match '^\s*\[(.+)\]\s*$') {
            if ($inWsl2 -and -not $done) { $out.Add("memory=${gb}GB"); $done = $true }
            $inWsl2 = ($Matches[1].Trim() -ieq 'wsl2'); if ($inWsl2) { $sawWsl2 = $true }
            $out.Add($l); continue
        }
        if ($inWsl2 -and $l -match '^\s*memory\s*=') {
            if (-not $done) { $out.Add("memory=${gb}GB"); $done = $true }
            continue
        }
        $out.Add($l)
    }
    if ($inWsl2 -and -not $done) { $out.Add("memory=${gb}GB"); $done = $true }
    if (-not $sawWsl2) { $out.Add('[wsl2]'); $out.Add("memory=${gb}GB") }
    if (Test-Path $f) { Copy-Item $f "$f.bak" -Force }
    Write-TextLF $f (($out -join "`r`n") + "`r`n")
}

function Restart-DockerDesktop {
    Info 'Restarting WSL and Docker Desktop to apply the new memory limit...'
    Get-Process -Name 'Docker Desktop', 'com.docker.backend', 'com.docker.build' -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
    Start-Sleep -Seconds 3
    & wsl.exe --shutdown
    Start-Sleep -Seconds 5
    Ensure-Docker
}

# Offered by setup: Docker Desktop's WSL 2 VM gets half the RAM by default,
# which on a 16 GB laptop caps the analysis at ~7 GB instead of 12 GB.
function Offer-DockerMemory {
    $h = Get-HostResources
    $d = Get-DockerResources
    $p = Get-PolicyBudget $h
    $targetMB = $p.memMB
    $useful = [int][math]::Round((Get-Requirement 'max_useful_memory_gb' 0) * 1024)
    if ($useful -gt 0) { $targetMB = [math]::Min($targetMB, $useful + 512) }   # + Docker's own share
    $targetGB = [int][math]::Ceiling($targetMB / 1024)
    Say ''
    Say "This computer: $($h.cpus) CPU cores, $(Format-GB $h.memMB) memory. Docker Desktop may use: $($d.cpus) cores, $(Format-GB $d.memMB)."
    if ($targetGB -lt 1 -or ($d.memMB + 1024) -ge ($targetGB * 1024)) { return }
    Say "The pipeline could use up to $targetGB GB (never more than 75% of the memory, leaving at least 3 GB to Windows), but Docker Desktop"
    Say "currently stops it at $(Format-GB $d.memMB). Raising Docker's limit writes memory=${targetGB}GB to $env:USERPROFILE\.wslconfig"
    Say '(a copy of any existing file is kept as .wslconfig.bak) and restarts Docker Desktop and WSL - anything else running'
    Say 'in WSL is stopped too. You can change it back at any time.'
    # --yes (agents, scripts) never answers this one: it restarts WSL, so it
    # needs the explicit --docker-memory, given after asking the user.
    $raise = $DockerMemory
    if (-not $raise -and -not $script:Yes) {
        try { $raise = ((Read-Host "Let Docker Desktop use up to $targetGB GB? [y/N]") -match '^(y|yes)$') } catch { $raise = $false }
    }
    if (-not $raise -and $script:Yes) {
        Say "Left unchanged. To raise it: run_crisprscreen setup --docker-memory"
        return
    }
    if ($raise) {
        Set-WslMemory $targetGB
        Restart-DockerDesktop
        $d = Get-DockerResources
        Say "Docker Desktop may now use $(Format-GB $d.memMB)."
    } else {
        Say 'Left unchanged. (Run setup again any time to change it.)'
    }
}

function Get-PipelineGit {
    $r = @{ commit = ''; modified = 'false' }
    if (-not (Have 'git')) { return $r }
    $c = & git -c safe.directory=* -C $PipelineDir rev-parse --short HEAD 2>$null
    if ($LASTEXITCODE -eq 0 -and $c) {
        $r.commit = "$c".Trim()
        $s = & git -c safe.directory=* -C $PipelineDir status --porcelain -- . 2>$null
        if ($s) { $r.modified = 'true' }
    }
    return $r
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
CRISPR_SCREEN pipeline (Docker version for Windows)

  run_crisprscreen setup            one-time: Docker Desktop + the pipeline's Docker image
  run_crisprscreen init [FOLDER]    set up a new project folder (default: current folder)
  run_crisprscreen                  run the pipeline for the project in the current folder
  run_crisprscreen run CONFIG.yaml  ... with another config file
  run_crisprscreen status           progress of the current/last run (and why it failed)
  run_crisprscreen logs             follow the live log
  run_crisprscreen stop             stop the running analysis (re-running resumes it)
  run_crisprscreen clean            delete work/ now (normally removed after a successful run unless keep_work: true)

Options: --detach  --yes  --pause  --rebuild  --docker-memory (setup)
Pipeline folder: $PipelineDir
"@
}

function Cmd-Setup {
    Info 'CRISPR_SCREEN setup (Docker version for Windows)'
    Ensure-Docker -AllowInstall
    Repair-LineEndings
    Ensure-Image -Force:$Rebuild

    # Docker Desktop (WSL 2) may use half the computer's memory by default: offer the budget's 75%.
    Offer-DockerMemory
    $p = Get-PolicyBudget (Get-HostResources)
    $d = Get-DockerResources
    $cpus = [math]::Min($p.cpus, $d.cpus); $mem = [math]::Min($p.memMB, $d.memMB - 512)
    $useful = [int][math]::Round((Get-Requirement 'max_useful_memory_gb' 0) * 1024)
    if ($useful -gt 0) { $mem = [math]::Min($mem, $useful) }
    Say "Each analysis will use at most $cpus CPU cores and $(Format-GB $mem) of memory (max_cpus / max_memory_gb in a project's config.yaml change this)."
    $minGB = Get-Requirement 'min_memory_gb' 2.5
    if ($mem -lt $minGB * 1024) {
        Warn "that is less than the $minGB GB this pipeline needs: runs will refuse to start."
    }

    # Tool versions, straight from the image.
    Say ''
    # No double quotes inside the bash command: Windows PowerShell 5.1 would strip them.
    & docker run --rm --platform linux/amd64 $Image bash -c 'echo cutadapt $(cutadapt --version); echo MAGeCK $(mageck -v 2>&1 | tail -1); fastqc --version; multiqc --version; python3 --version; R --version | head -1; nextflow -version 2>&1 | grep -i version'

    # Put run_crisprscreen on the user's PATH (new terminals only).
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if (-not $userPath) { $userPath = '' }
    $parts = $userPath -split ';' | Where-Object { $_ }
    if (-not ($parts | Where-Object { $_.TrimEnd('\') -ieq $PipelineDir.TrimEnd('\') })) {
        [Environment]::SetEnvironmentVariable('Path', (($parts + $PipelineDir) -join ';'), 'User')
        Say ''
        Say "Added $PipelineDir to your PATH: in a NEW terminal window, 'run_crisprscreen' works from any folder."
    }
    Say ''
    Info 'Setup complete. Next: create a project folder with   run_crisprscreen init <folder>'
    Done 0
}

function Cmd-Init {
    $target = (Get-Location).ProviderPath
    if ($positional.Count -gt 0) { $target = $positional[0] }
    New-Item -ItemType Directory -Force $target | Out-Null
    $target = (Resolve-Path -LiteralPath $target).ProviderPath.TrimEnd('\')
    $repoRoot = Split-Path $PipelineDir -Parent
    if ($target.StartsWith($repoRoot, [StringComparison]::OrdinalIgnoreCase)) {
        Fail "Don't put a project inside the pipeline's own folder ($repoRoot). Choose another folder, e.g. $HOME\Documents\my_screen"
    }
    Info "Setting up a new $PipelineName project in $target ..."
    foreach ($d in @('original_data\fastq', 'metadata', 'results', 'utils')) {
        New-Item -ItemType Directory -Force (Join-Path $target $d) | Out-Null
    }
    $cfg = Join-Path $target 'config.yaml'
    if (Test-Path $cfg) { Say '  config.yaml already exists - leaving it alone.' }
    else { Copy-Item (Join-Path $PipelineDir 'examples\config.yaml') $cfg; Say '  Fetched config.yaml' }

    $ex = Join-Path $target 'examples'
    New-Item -ItemType Directory -Force $ex | Out-Null
    foreach ($f in Get-ChildItem (Join-Path $PipelineDir 'examples') -File) {
        $d = Join-Path $ex $f.Name
        if (-not (Test-Path $d)) { Copy-Item $f.FullName $d }
    }
    Say '  Fetched examples\ (samplesheet_example.csv, library_example.csv, contrasts_example.csv)'

    # Default barcode-position plot (utils\): fetched once, never overwritten. Its md5 goes
    # into utils\.pipeline_defaults.md5 (md5sum format) so later changes can be spotted.
    $srcPlot = Join-Path $PipelineDir 'bin\plot_barcode_position.R'
    $dstPlot = Join-Path $target 'utils\plot_barcode_position.R'
    if (Test-Path -LiteralPath $srcPlot) {
        if (Test-Path -LiteralPath $dstPlot) {
            Say '  utils\plot_barcode_position.R already exists - leaving it alone.'
        } else {
            Copy-Item -LiteralPath $srcPlot -Destination $dstPlot
            $md5 = (Get-FileHash -Algorithm MD5 -LiteralPath $dstPlot).Hash.ToLower()
            $md5File = Join-Path $target 'utils\.pipeline_defaults.md5'
            $kept = @()
            if (Test-Path -LiteralPath $md5File) {
                $kept = @([System.IO.File]::ReadAllLines($md5File) | Where-Object { -not $_.EndsWith('  utils/plot_barcode_position.R') })
            }
            $kept += "$md5  utils/plot_barcode_position.R"
            Write-TextLF $md5File (($kept -join "`n") + "`n")
            Say '  Fetched utils\plot_barcode_position.R (default barcode-position plot; md5 recorded in utils\.pipeline_defaults.md5)'
        }
    }

    # AI agent settings (.claude\) from the pipeline's agent\ folder; existing files are never overwritten.
    $srcClaude = Join-Path $PipelineDir 'agent\.claude'
    if (Test-Path -LiteralPath $srcClaude -PathType Container) {
        $dstClaude = Join-Path $target '.claude'
        $srcRoot = $srcClaude.TrimEnd('\') + '\'
        $copied = 0
        foreach ($f in Get-ChildItem -LiteralPath $srcClaude -Recurse -File -Force) {
            $rel = $f.FullName.Substring($srcRoot.Length)
            $dst = Join-Path $dstClaude $rel
            if (Test-Path -LiteralPath $dst) { continue }
            New-Item -ItemType Directory -Force (Split-Path $dst -Parent) | Out-Null
            Copy-Item -LiteralPath $f.FullName -Destination $dst
            $copied++
        }
        if ($copied -gt 0) { Say "  Fetched .claude\ ($copied files for AI agents; existing files left alone)" }
    }

    # Agent instructions, with this computer's pipeline location filled in.
    foreach ($f in @('AGENTS.md', 'CLAUDE.md')) {
        $d = Join-Path $target $f
        if (Test-Path $d) { Say "  $f already exists - leaving it alone."; continue }
        $text = [System.IO.File]::ReadAllText((Join-Path $PipelineDir "agent\$f"))
        Write-TextLF $d ($text.Replace('{{PIPELINE_DIR}}', $PipelineDir))
        Say "  Fetched $f"
    }

    # Launchers inside the project (rewritten each time: they hold the pipeline's location).
    # The terminal launcher lives in utils\ (the project root keeps only the double-click RUN_ANALYSIS.bat).
    $bat = Join-Path $PipelineDir 'run_crisprscreen.bat'
    $marker = 'rem Runs the CRISPR_SCREEN pipeline (Docker version) for THIS project folder.'
    New-Item -ItemType Directory -Force (Join-Path $target 'utils') | Out-Null
    $cmd = "@echo off`r`n$marker`r`n" +
           "rem Same commands as the pipeline's own launcher: (no argument) = run, status, logs, stop, clean.`r`n" +
           "rem Generated by 'run_crisprscreen init'; the pipeline lives in:`r`nrem   $PipelineDir`r`n" +
           "cd /d `"%~dp0..`"`r`ncall `"$bat`" %*`r`nexit /b %ERRORLEVEL%`r`n"
    [System.IO.File]::WriteAllText((Join-Path $target 'utils\run_crisprscreen.bat'), $cmd, [System.Text.Encoding]::ASCII)
    # An older init wrote the terminal launcher to the project root: remove it only if it is ours (marker line).
    $oldBat = Join-Path $target 'run_crisprscreen.bat'
    if ((Test-Path -LiteralPath $oldBat -PathType Leaf) -and ([System.IO.File]::ReadAllText($oldBat).Contains($marker))) {
        Remove-Item -LiteralPath $oldBat -Force
        Say 'Removed the old run_crisprscreen.bat from the project root (now utils\run_crisprscreen.bat)'
    }
    $dbl = "@echo off`r`nrem Double-click to run (or resume) the CRISPR screen analysis of this folder.`r`n" +
           "cd /d `"%~dp0`"`r`ncall `"$bat`" run --pause`r`n"
    [System.IO.File]::WriteAllText((Join-Path $target 'RUN_ANALYSIS.bat'), $dbl, [System.Text.Encoding]::ASCII)
    Say '  Created utils\run_crisprscreen.bat (terminal) and RUN_ANALYSIS.bat (double-click)'

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
        $null = & git -C $target -c user.name=crisprscreen -c user.email=crisprscreen@localhost commit -q -m "Initialize $PipelineName project" 2>&1
        Say '  Initialized git repository (raw data / results / work excluded - see .gitignore; utils\ and .claude\ are included in the first commit)'
    } else {
        Say '  (git not installed - skipped version control; optional)'
    }

    Say ''
    Say 'Next steps:'
    Say '  1. Put your input tables in this folder (start from the files in examples\):'
    Say '       samplesheet.csv, library.csv, and - when you are ready for mageck test - contrasts.csv'
    Say '  2. Edit config.yaml: project_name, and fastq_r1 (multiplexed FASTQ) or fastq_dir (folder of already-demultiplexed FASTQs).'
    Say '     Your FASTQ files stay where they are: put their full path in config.yaml and keep it stable.'
    Say '  3. Run:  double-click RUN_ANALYSIS.bat   (or type  .\utils\run_crisprscreen.bat  in this folder)'
    Done 0
}

function Cmd-Run {
    Get-ProjectDirFromArgs
    if (-not (Test-Path -LiteralPath $script:ConfigFile)) {
        Fail "$($script:ConfigName) not found in $($script:ProjectDir). Set the project up first:  run_crisprscreen init"
    }
    $name = Get-ContainerName $script:ProjectDir
    Ensure-Docker
    if (Test-Running $name) {
        Fail "this project's analysis is already running. Watch it with:  run_crisprscreen logs   (or stop it: run_crisprscreen stop)"
    }
    Repair-LineEndings
    $budget = Get-RunBudget
    if ((Get-DriveMediaType $script:ProjectDir) -eq 'HDD') {
        Warn ("this project folder is on a hard disk drive (HDD). The analysis works, but demultiplexing will be several " +
              "times slower than on an SSD, and Windows may respond slowly meanwhile. An internal SSD is better if it has room.")
    }

    Info "Checking the input files named in $($script:ProjectDir)\$($script:ConfigName)"
    # FASTQs outside the project folder: read-only mounts, given to the container below.
    $inputArgs = @()
    $fqDir = Resolve-RawFastq 'fastq_dir'
    if ($null -ne $fqDir) {
        $r1 = Get-YamlValue $script:ConfigFile 'fastq_r1'
        if (-not (Test-Unset $r1)) { Say '  fastq_dir is set, so the samples are treated as already demultiplexed (fastq_r1 is ignored).' }
        $inputArgs += @(Get-FastqInputArgs $fqDir)
    } else {
        $r1 = Get-YamlValue $script:ConfigFile 'fastq_r1'
        if (Test-Unset $r1) { Fail 'config.yaml needs fastq_r1 (a multiplexed FASTQ) or fastq_dir (a folder of already-demultiplexed FASTQs).' }
        $inputArgs += @(Get-FastqInputArgs (Resolve-RawFastq 'fastq_r1'))
        $fqR2 = Resolve-RawFastq 'fastq_r2'
        if ($null -ne $fqR2) { $inputArgs += @(Get-FastqInputArgs $fqR2) }
    }
    Check-ConfigPaths

    Ensure-Image

    $ts = Get-Date -Format 'yyyyMMdd_HHmmss'
    $git = Get-PipelineGit
    $imageId = "$(& docker image inspect $Image --format '{{.Id}}' 2>$null)".Trim()
    $dockerVersion = "$(& docker version --format '{{.Server.Version}}' 2>$null)".Trim()
    $os = 'Windows'
    try { $o = Get-CimInstance Win32_OperatingSystem; $os = "$($o.Caption) $($o.Version)".Trim() } catch { }
    $configRel = $script:ConfigName
    # Mark the run as started now: until the container rewrites it, `status`
    # (and an agent polling it) would otherwise still see the previous run's outcome.
    $logDir = Join-Path $script:ProjectDir 'metadata\nf_logs'
    New-Item -ItemType Directory -Force $logDir | Out-Null
    Write-TextLF (Join-Path $logDir 'last_run.status') "running`n"

    $dockerArgs = @(
        'run', '-d', '--rm', '--init', '--name', $name, '--platform', 'linux/amd64',
        '--label', "crisprscreen.project=$($script:ProjectDir)",
        '--memory', "$($budget.memMB)m", '--memory-swap', "$($budget.memMB)m",
        '--cpuset-cpus', "0-$($budget.cpus - 1)",
        '-e', "CRISPR_BUDGET_CPUS=$($budget.cpus)",
        '-e', "CRISPR_BUDGET_MEMORY_MB=$($budget.memMB)",
        '-e', "CRISPR_BUDGET_REASON=$($budget.reason)",
        '--mount', "type=bind,source=$($script:ProjectDir),target=/project",
        '--mount', "type=bind,source=$PipelineDir,target=/pipeline,readonly"
    )
    $dockerArgs += $inputArgs
    $dockerArgs += @(
        '-e', "CRISPR_IMAGE=$Image",
        '-e', "CRISPR_IMAGE_ID=$imageId",
        '-e', 'CRISPR_CONTAINER_ENGINE=docker',
        '-e', "CRISPR_CONTAINER_ENGINE_VERSION=$dockerVersion",
        '-e', "CRISPR_HOST_OS=Docker Desktop on $os",
        '-e', "CRISPR_PIPELINE_COMMIT=$($git.commit)",
        '-e', "CRISPR_PIPELINE_MODIFIED=$($git.modified)",
        '-e', "CRISPR_PIPELINE_HOST_DIR=$PipelineDir",
        '-e', "CRISPR_PROJECT_HOST_DIR=$($script:ProjectDir)",
        $Image, 'bash', '/pipeline/docker/run_in_container.sh', $configRel, $ts
    )
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
    Say  '  Progress : run_crisprscreen status      Stop: run_crisprscreen stop'
    Say  '  Closing this window does NOT stop the analysis. Keep the computer awake until it finishes.'
    if ($Detach) { Done 0 }
    Say  '  (Press Ctrl+C to stop watching; the analysis keeps running.)'
    Say ''
    Follow-Run $name
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

# Hidden helper for as long as a run's container (or, during an image build,
# this launcher's process) lives:
#  - keeps Windows awake (no sleep mid-analysis);
#  - guards Windows' memory: Docker's Linux VM keeps the files it has read
#    (FASTQs, image layers) as cache, up to everything Docker may use, and
#    hands that memory back to Windows only slowly. On a computer whose memory
#    is already busy that pushes Windows into paging until it stops responding.
#    When Windows' free memory drops below max(4 GB, 20% of RAM), the helper
#    tells the VM to drop that cache (clean, re-readable data; nothing is lost).
function Start-Guard([string]$container = '', [int]$watchPid = 0) {
    $totalMB = 0
    try { $totalMB = [int]((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1MB) } catch { }
    $minFreeMB = [math]::Max(4096, [int]($totalMB / 5))
    $alive = "[bool](Get-Process -Id $watchPid -ErrorAction SilentlyContinue)"
    if ($container) { $alive = "[bool](& docker ps -q --filter `"name=^/$container`$`" 2>`$null)" }
    $code = @"
Add-Type -Namespace CrisprScreen -Name Power -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint flags);'
`$i = 0
while ($alive) {
    if (`$i % 6 -eq 0) { [CrisprScreen.Power]::SetThreadExecutionState([uint32]2147483649) | Out-Null }
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
        Say  "  Results : $($script:ProjectDir)\results\   (start with results\counts\qc\report\*_qc_report.html)"
        Say  "  Methods : $($script:ProjectDir)\METHODS.md"
        Done 0
    } elseif ($st -eq 'running' -and (Test-Running $name)) {
        Say 'Stopped watching; the analysis is still running.'
        Done 0
    } else {
        Write-Host "Pipeline did not finish: $st" -ForegroundColor Red
        Say '  See:  run_crisprscreen status   (it names the failed step and where its messages were saved)'
        Say '  Fix the problem, then run again: finished steps are reused, not redone.'
        Done 1
    }
}

function Cmd-Status {
    Get-ProjectDirFromArgs
    $name = Get-ContainerName $script:ProjectDir
    $running = $false
    if ((Have 'docker') -and ((Quiet { docker info }) -eq 0)) { $running = Test-Running $name }
    $statusFile = Join-Path $script:ProjectDir 'metadata\nf_logs\last_run.status'
    $logs = @(Get-ChildItem (Join-Path $script:ProjectDir 'metadata\nf_logs') -Filter 'nextflow_*.log' -ErrorAction SilentlyContinue | Sort-Object Name)
    if ($logs.Count -eq 0) { Say "No run yet in $($script:ProjectDir)."; Done 0 }
    $last = $logs[-1]
    if ($running) { Info 'Status: RUNNING' }
    elseif (Test-Path $statusFile) {
        $st = ([System.IO.File]::ReadAllText($statusFile)).Trim()
        if ($st -eq 'running') { $st = 'interrupted (Docker or the computer stopped mid-run; run again to resume)' }
        Info "Status: last run $st"
    }
    Say "Log: $($last.FullName)"

    $trace = Join-Path $script:ProjectDir 'metadata\pipeline_info\execution_trace.txt'
    if (Test-Path $trace) {
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
    $work = Join-Path $script:ProjectDir 'work'
    $nxDir = Join-Path $script:ProjectDir '.nextflow'
    if (-not (Test-Path $work) -and -not (Test-Path $nxDir)) { Say 'No work\ folder - nothing to clean.'; Done 0 }
    Ensure-Docker
    if (Test-Running $name) { Fail 'the analysis is running - stop it first (run_crisprscreen stop).' }
    Say "This deletes $work : Nextflow's working copies of every step."
    Say 'results\, original_data\ and METHODS.md are NOT touched, but the next run will start from scratch.'
    if (-not (Confirm 'Delete it?')) { Say 'Nothing deleted.'; Done 0 }
    # Deleted from inside Linux: work\ holds Linux symlinks that Windows tools handle poorly.
    $null = & docker run --rm --platform linux/amd64 --mount "type=bind,source=$($script:ProjectDir),target=/project" $Image rm -rf /project/work /project/.nextflow 2>&1
    if (Test-Path $work) { Remove-Item -LiteralPath $work -Recurse -Force -ErrorAction SilentlyContinue }
    Say 'Deleted.'
    Done 0
}

switch -Regex ($Command) {
    '^setup$'  { Cmd-Setup }
    '^init$'   { Cmd-Init }
    '^run$'    { Cmd-Run }
    '^status$' { Cmd-Status }
    '^logs$'   { Cmd-Logs }
    '^stop$'   { Cmd-Stop }
    '^clean$'  { Cmd-Clean }
    default    { Show-Help; Done 0 }
}
