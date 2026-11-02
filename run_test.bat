@echo off
rem Checks that Python can drive DWSIM. Expected last line: ... PASS
call "%~dp0run.bat" test_dwsim.py
