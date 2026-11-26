@echo off
rem famscan launcher for Windows: no venv activation and no "python" command needed.
rem From this folder, in PowerShell:  .\famscan start    (also: .\famscan demo, .\famscan list, ...)
rem In Command Prompt the .\ is optional:  famscan start
if exist "%~dp0.venv\Scripts\python.exe" (
  "%~dp0.venv\Scripts\python.exe" -m famscan.cli %*
  exit /b %errorlevel%
)
echo famscan is not set up in this folder yet. Run these two lines once:
echo     py -3.12 -m venv .venv
echo     .venv\Scripts\python -m pip install -e .
echo (No "py"? Install Python 3.12 from python.org, ticking "Add python.exe to PATH".)
echo Then:  .\famscan start
exit /b 1
