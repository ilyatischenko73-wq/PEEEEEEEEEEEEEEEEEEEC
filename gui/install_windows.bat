@echo off
setlocal
pushd "%~dp0.."
py -3 -m venv .venv
if errorlevel 1 goto failed
.venv\Scripts\python.exe -m pip install -r gui\requirements.txt
if errorlevel 1 goto failed
echo Installation complete. Run gui\start_windows.bat
popd
pause
exit /b 0
:failed
echo Installation failed. Check Python installation and the messages above.
popd
pause
exit /b 1
