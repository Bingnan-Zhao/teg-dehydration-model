@echo off
rem Finds the Anaconda/Miniconda install and sets CONDA_ROOT.
set "CONDA_ROOT="
for %%D in ("%USERPROFILE%\anaconda3" "%USERPROFILE%\miniconda3" "%LOCALAPPDATA%\anaconda3" "%LOCALAPPDATA%\miniconda3" "%ProgramData%\anaconda3" "%ProgramData%\miniconda3" "C:\anaconda3" "C:\miniconda3") do (
  if not defined CONDA_ROOT if exist "%%~D\Scripts\activate.bat" set "CONDA_ROOT=%%~D"
)
if not defined CONDA_ROOT (
  for /f "delims=" %%P in ('where conda 2^>nul') do if not defined CONDA_ROOT (
    for %%Q in ("%%~dpP..") do if exist "%%~fQ\Scripts\activate.bat" set "CONDA_ROOT=%%~fQ"
  )
)
if not defined CONDA_ROOT (
  echo [ERROR] Anaconda not found. Edit _find_conda.bat and set CONDA_ROOT to your Anaconda folder.
  exit /b 1
)
call "%CONDA_ROOT%\Scripts\activate.bat" "%CONDA_ROOT%"
exit /b 0
