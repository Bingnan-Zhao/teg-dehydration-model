@echo off
rem One-time setup: creates the conda environment "teg" with the packages this project needs.
cd /d "%~dp0"
call _find_conda.bat || (pause & exit /b 1)
echo Using Anaconda at %CONDA_ROOT%
call conda env list | findstr /b /c:"teg " >nul
if errorlevel 1 (
  echo Creating environment "teg" ^(Python 3.12^)...
  call conda create -y -n teg python=3.12 || (echo [ERROR] conda create failed & pause & exit /b 1)
) else (
  echo Environment "teg" already exists, updating packages...
)
call conda activate teg
python -m pip install --upgrade pip
python -m pip install pythonnet numpy pandas scipy matplotlib fmpy || (echo [ERROR] pip install failed & pause & exit /b 1)
echo.
echo Setup finished. Next: double-click run_test.bat
pause
