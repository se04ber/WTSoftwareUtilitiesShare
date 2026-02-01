@echo off
setlocal enableextensions

REM Start Jupyter Notebook from the dedicated WTConc2 venv.
REM This avoids launching Anaconda's Jupyter by accident.

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..\..\..\") do set "REPO_ROOT=%%~fI"
set "PROJECT_DIR=%REPO_ROOT%\Manual\WTSoftwareUtilitiesShare"
set "ENV_DIR=%PROJECT_DIR%\.venv_WTConc2"
set "VENV_PY=%ENV_DIR%\Scripts\python.exe"

if not exist "%VENV_PY%" (
  echo ERROR: WTConc2 venv not found. Run setup first:
  echo "%SCRIPT_DIR%setup_WTConc2_venv.bat"
  exit /b 1
)

pushd "%PROJECT_DIR%"

REM Ensure local imports work (windtunnel package lives in this folder)
set PYTHONNOUSERSITE=1

echo Starting Jupyter in: "%CD%"
echo Using python: "%VENV_PY%"
echo.

call "%VENV_PY%" -m notebook

popd


