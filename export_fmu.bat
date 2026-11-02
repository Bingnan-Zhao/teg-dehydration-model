@echo off
rem Step 8 (Dymola): simulate the step test and export modelica\TEGUnit.fmu (FMI 2.0, ME + CS)
setlocal
set "DYMOLA="
for /d %%D in ("%ProgramFiles%\Dymola*") do if exist "%%~D\bin64\Dymola.exe" set "DYMOLA=%%~D\bin64\Dymola.exe"
if not defined DYMOLA (
  echo [ERROR] Dymola not found under "%ProgramFiles%". Edit export_fmu.bat and set DYMOLA to your Dymola.exe.
  pause
  exit /b 1
)
echo Using %DYMOLA%
set "P=%~dp0modelica"
set "P=%P:\=/%"
(
  echo cd^("%P%"^);
  echo openModel^("%P%/TEGDehydration/package.mo"^);
  echo cd^("%P%"^);
  echo simulateModel^("TEGDehydration.Examples.StepTest", stopTime=93600, outputInterval=60, method="Dassl", tolerance=1e-6, resultFile="StepTest"^);
  echo translateModelFMU^("TEGDehydration.TEGUnit", false, "TEGUnit", "2", "all", false^);
  echo exit^(^);
) > "%~dp0modelica\_export.mos"
"%DYMOLA%" /nowindow "%~dp0modelica\_export.mos"
echo.
if exist "%~dp0modelica\TEGUnit.fmu" (echo FMU exported: modelica\TEGUnit.fmu) else (echo [WARNING] TEGUnit.fmu not found. Open modelica\TEGDehydration\package.mo in Dymola and check the messages.)
pause
