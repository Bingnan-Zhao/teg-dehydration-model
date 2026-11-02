@echo off
rem Step 3: fit the Peng-Robinson TEG-water kij to open VLE data (about 1-2 minutes)
call "%~dp0run.bat" 03_fit_teg_water.py
