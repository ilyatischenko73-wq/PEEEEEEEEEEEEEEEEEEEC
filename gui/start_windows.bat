@echo off
setlocal
pushd "%~dp0.."
if not exist .venv\Scripts\python.exe (
  echo First run gui\install_windows.bat
  popd
  pause
  exit /b 1
)
.venv\Scripts\python.exe peec_gui.py
set "peec_gui_result=%errorlevel%"
if not "%peec_gui_result%"=="0" pause
popd
exit /b %peec_gui_result%
