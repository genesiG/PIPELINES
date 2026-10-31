@echo off
rem run_cutandrun.bat - CUT&RUN pipeline launcher (Docker version for Windows).
rem All the work is in run_cutandrun.ps1 next to this file; this wrapper only
rem starts it with an execution policy that lets it run on any Windows computer.
rem   run_cutandrun setup | init [FOLDER] | [run] [CONFIG] | status | logs | stop | clean | script utils\X.R [ARGS] | help
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_cutandrun.ps1" %*
exit /b %ERRORLEVEL%
