@echo off
setlocal
where ifx >nul 2>nul
if errorlevel 1 (
  echo Open Intel oneAPI command prompt for Visual Studio, then run this script.
  exit /b 1
)
pushd "%~dp0..\.."
if not exist "build\windows" mkdir "build\windows"
if not exist "build\windows" exit /b 1
pushd "build\windows"
ifx /nologo /O2 /Qopenmp /Qmkl:sequential /module:. ^
  "..\..\fortran\src\peec_base.f90" ^
  "..\..\fortran\src\peec_linalg.f90" ^
  "..\..\fortran\src\peec_quadrature.f90" ^
  "..\..\fortran\src\peec_mesh.f90" ^
  "..\..\fortran\src\peec_config.f90" ^
  "..\..\fortran\src\peec_model.f90" ^
  "..\..\fortran\src\peec_physics.f90" ^
  "..\..\fortran\src\peec_transient.f90" ^
  "..\..\fortran\src\peec_aperture_json.f90" ^
  "..\..\fortran\src\peec_io.f90" ^
  "..\..\fortran\app\main.f90" /exe:PeecSolverFortran.exe
set "peec_result=%errorlevel%"
popd
popd
exit /b %peec_result%
