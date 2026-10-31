@echo off
rem run_rnaseq.bat - RNA-seq pipeline launcher (Docker version for Windows).
rem All the work is in run_rnaseq.ps1 next to this file; this wrapper only
rem starts it with an execution policy that lets it run on any Windows computer.
rem   run_rnaseq setup | init [FOLDER] | [run] [CONFIG] | status | logs | stop | clean | script utils\featurecounts.sh [ARGS] | help
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_rnaseq.ps1" %*
exit /b %ERRORLEVEL%
