@echo off
rem ============================================================
rem  PDOS Structure Viewer - GUI launcher
rem  Left: projects | Middle: structure/DOS | Right: parameters
rem ============================================================
setlocal EnableExtensions
cd /d "%~dp0"

set "SCRIPT=%~dp0pdos_viewer.py"
set "PY=python.exe"

if not exist "%SCRIPT%" (
    echo [ERROR] pdos_viewer.py not found:
    echo         %SCRIPT%
    pause
    exit /b 1
)

for %%I in (python.exe) do set "PY=%%~$PATH:I"
if not exist "%PY%" set "PY=python.exe"

rem prefer pythonw.exe (no console); fall back to python.exe
set "PYW=%PY:python.exe=pythonw.exe%"
if not exist "%PYW%" set "PYW=%PY%"

start "" "%PYW%" "%SCRIPT%"
exit /b 0
