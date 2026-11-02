@echo off
rem Step 2: tune the Peng-Robinson methane-water kij (1-2 minutes)
call "%~dp0run.bat" 02_tune_kij.py
