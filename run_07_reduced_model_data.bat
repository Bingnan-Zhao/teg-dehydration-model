@echo off
rem Step 7: fit the reduced-order model relations from DWSIM (3-5 minutes; coefficients are already included)
call "%~dp0run.bat" 07_reduced_model_data.py
