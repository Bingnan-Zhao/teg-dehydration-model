@echo off
rem Runs a script from the scripts folder in the "teg" environment.
rem Usage: run.bat script_name.py   (or drag a .py file onto this .bat)
cd /d "%~dp0"
if "%~1"=="" (echo Usage: run.bat script_name.py & pause & exit /b 1)
call _find_conda.bat || (pause & exit /b 1)
call conda activate teg || (echo [ERROR] Run setup_env.bat first & pause & exit /b 1)
cd /d "%~dp0scripts"
python "%~nx1"
echo.
pause
