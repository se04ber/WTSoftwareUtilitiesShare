
REM Minimal Jupyter launcher for WTConc2

cd /d "%~dp0"
set "JUP=%CD%\WTConc2\Scripts\jupyter.exe"

IF NOT EXIST "%JUP%" (
    echo WTConc2 not found. Run windows_deploy\setup_WTConc2.bat first.
)

"%JUP%" notebook
