
REM Minimal setup script for WTConc2 (Batch)
REM Prefer Python 3.11, else default



set "PY_CMD=C:\Users\u302117\AppData\Local\Programs\Python\Launcher\py.exe"
set "ENV_DIR=%CD%\WTConc3"
set "VENV_PY=%ENV_DIR%\Scripts\python.exe"


%PY_CMD% -m venv "%ENV_DIR%"

IF NOT EXIST "%VENV_PY%" (
    echo WTConc3 venv python not found.
    
)

IF NOT EXIST "%VENV_PY%" (
    echo WTConc3 venv python not found.
    exit /b 1
) ELSE (
    call "%ENV_DIR%\Scripts\activate.bat"
)

set PYTHONNOUSERSITE=1
set PIP_REQUIRE_VIRTUALENV=1

echo "%VENV_PY%"
echo "%CD%"

pip install --upgrade pip setuptools wheel

IF EXIST "%CD%\requirements.txt" (
    pip install -r "%CD%\requirements.txt"
)
pip install "git+https://github.com/se04ber/WTSoftwareUtilitiesShare.git"
pip install notebook ipykernel
ipykernel install --user --name WTConc3 --display-name "Python (WTConc3)"

echo WTConc2 ready. Activate with: "%ENV_DIR%\Scripts\activate.bat"
echo Start Jupyter with: "%ENV_DIR%\Scripts\jupyter.exe" notebook

cmd /k