@echo off
rem run_atacseq.bat - ATAC-seq pipeline launcher (Docker version for Windows).
rem All the work is in run_atacseq.ps1 next to this file; this wrapper only
rem starts it with an execution policy that lets it run on any Windows computer.
rem   run_atacseq setup | init [FOLDER] | [run] [CONFIG] | status | logs | stop | clean | script utils\X.R [ARGS] | help
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_atacseq.ps1" %*
exit /b %ERRORLEVEL%
