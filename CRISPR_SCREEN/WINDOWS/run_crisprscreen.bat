@echo off
rem run_crisprscreen.bat - CRISPR_SCREEN pipeline launcher (Docker version for Windows).
rem All the work is in run_crisprscreen.ps1 next to this file; this wrapper only
rem starts it with an execution policy that lets it run on any Windows computer.
rem   run_crisprscreen setup | init [FOLDER] | [run] [CONFIG] | status | logs | stop | clean | script utils\X.R [ARGS] | help
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_crisprscreen.ps1" %*
exit /b %ERRORLEVEL%
