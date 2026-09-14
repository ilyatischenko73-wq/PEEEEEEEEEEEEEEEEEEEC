@echo off
setlocal
pushd "%~dp0..\.."
if not exist "build\windows\PeecSolverFortran.exe" (
  echo First run fortran\windows\build_ifx.bat from Intel oneAPI command prompt.
  popd
  exit /b 1
)
"build\windows\PeecSolverFortran.exe" %*
set "peec_result=%errorlevel%"
pause
popd
exit /b %peec_result%
