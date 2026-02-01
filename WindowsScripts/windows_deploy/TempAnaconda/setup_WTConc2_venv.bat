@echo off
setlocal enableextensions enabledelayedexpansion

REM =============================================================================
REM WTConc2 venv setup for Windows (works alongside Anaconda).
REM
REM Why: avoid PATH/kernel/package mixing between Anaconda and this project.
REM
REM Usage:
REM   setup_WTConc2_venv.bat
REM     - creates/updates venv, installs deps, registers kernel "WTConc2 (venv)"
REM
REM   setup_WTConc2_venv.bat download-wheels
REM     - downloads wheels into ".\wheels\" for offline install (run on online PC)
REM =============================================================================

REM Script directory (TempAnaconda)
set "SCRIPT_DIR=%~dp0"

REM Repo root: TempAnaconda -> windows_deploy -> WindowsScripts -> repo root
for %%I in ("%SCRIPT_DIR%..\..\..\") do set "REPO_ROOT=%%~fI"

set "PROJECT_DIR=%REPO_ROOT%\Manual\WTSoftwareUtilitiesShare"
set "REQ_MAIN=%PROJECT_DIR%\requirements.txt"
set "REQ_DEV=%PROJECT_DIR%\requirements-dev.txt"

set "ENV_DIR=%PROJECT_DIR%\.venv_WTConc2"
set "VENV_PY=%ENV_DIR%\Scripts\python.exe"

echo Repo root:     "%REPO_ROOT%"
echo Project dir:   "%PROJECT_DIR%"
echo Venv dir:      "%ENV_DIR%"
echo.

if not exist "%PROJECT_DIR%" (
  echo ERROR: Project folder not found: "%PROJECT_DIR%"
  echo This script must live inside the WTSoftwareUtilitiesShare repo.
  exit /b 1
)

REM Prefer python.org launcher (py). This avoids accidentally using conda python.
where py >nul 2>nul
if errorlevel 1 (
  echo ERROR: Windows "py" launcher not found.
  echo Please install Python from python.org (recommended: 3.11.x) and include the py launcher.
  exit /b 1
)

REM Pick a Python version.
REM We REQUIRE Python >= 3.11 for this setup (tested in CI); prefer 3.11, allow 3.12.
set "PY_CMD="

py -3.11 -c "import sys; print(sys.version)" >nul 2>nul
if not errorlevel 1 (
  set "PY_CMD=py -3.11"
) else (
  py -3.12 -c "import sys; print(sys.version)" >nul 2>nul
  if not errorlevel 1 (
    set "PY_CMD=py -3.12"
  )
)

if "%PY_CMD%"=="" (
  echo ERROR: Supported Python version not found.
  echo This setup requires Python 3.11+ (recommended 3.11.x).
  echo.
  echo Installed Python interpreters detected by the py launcher:
  py -0p
  echo.
  echo Fix: install Python 3.11 from python.org (and include the py launcher).
  exit /b 1
)

echo Using: %PY_CMD%
%PY_CMD% -c "import sys, platform; print('Python:', sys.version); print('Executable:', sys.executable); print('Arch:', platform.architecture()[0])"
echo.

REM Optional: download wheels for offline installation
if /I "%~1"=="download-wheels" goto :download_wheels

REM Create venv
%PY_CMD% -m venv "%ENV_DIR%"
if not exist "%VENV_PY%" (
  echo ERROR: venv python not found: "%VENV_PY%"
  exit /b 1
)

set PYTHONNOUSERSITE=1
set PIP_REQUIRE_VIRTUALENV=1

REM Install deps (offline-capable if wheels folder exists)
set "WHEEL_DIR=%SCRIPT_DIR%wheels"

call "%VENV_PY%" -m pip install --upgrade pip setuptools wheel

if exist "%WHEEL_DIR%" (
  echo Installing from local wheels: "%WHEEL_DIR%"
  if exist "%REQ_MAIN%" (
    call "%VENV_PY%" -m pip install --no-index --find-links "%WHEEL_DIR%" -r "%REQ_MAIN%"
  )
  if exist "%REQ_DEV%" (
    call "%VENV_PY%" -m pip install --no-index --find-links "%WHEEL_DIR%" -r "%REQ_DEV%"
  )
  call "%VENV_PY%" -m pip install --no-index --find-links "%WHEEL_DIR%" notebook ipykernel
) else (
  echo Installing from PyPI (internet required).
  echo If you need offline install: run this script with "download-wheels" on an online machine
  echo and copy the resulting "%WHEEL_DIR%" folder to the offline machine.
  if exist "%REQ_MAIN%" (
    call "%VENV_PY%" -m pip install -r "%REQ_MAIN%"
  )
  if exist "%REQ_DEV%" (
    call "%VENV_PY%" -m pip install -r "%REQ_DEV%"
  )
  call "%VENV_PY%" -m pip install notebook ipykernel
)

REM Register kernel for ANY Jupyter (including Anaconda's) to use
call "%VENV_PY%" -m ipykernel install --user --name WTConc2-venv --display-name "WTConc2 (venv)"

echo.
echo Done.
echo - To start Jupyter (recommended): run "%SCRIPT_DIR%start_WTConc2_Jupyter.bat"
echo - Or open Anaconda Jupyter and select kernel: "WTConc2 (venv)"
echo.
pause

goto :eof

:download_wheels
set "WHEEL_DIR=%SCRIPT_DIR%wheels"
echo Downloading wheels to: "%WHEEL_DIR%"
echo.
if not exist "%WHEEL_DIR%" mkdir "%WHEEL_DIR%"

REM Use the selected Python to download wheels for that interpreter.
%PY_CMD% -m pip download --dest "%WHEEL_DIR%" -r "%REQ_MAIN%"
if exist "%REQ_DEV%" (
  %PY_CMD% -m pip download --dest "%WHEEL_DIR%" -r "%REQ_DEV%"
)
%PY_CMD% -m pip download --dest "%WHEEL_DIR%" notebook ipykernel

echo.
echo Done. Copy "%WHEEL_DIR%" to the offline machine (same relative path).
echo Then run setup_WTConc2_venv.bat (without arguments) on the offline machine.
echo.
pause

