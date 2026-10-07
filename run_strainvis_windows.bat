@echo off
setlocal enabledelayedexpansion

:: Change working directory to the directory where this batch file lives
cd /d "%~dp0"

echo Searching for Conda installation...

:: 1. Check if CONDA_EXE is already set (e.g., launched inside Anaconda Prompt)
if defined CONDA_EXE (
    set "CONDA_BAT=%CONDA_EXE:\Scripts\conda.exe=\condabin\conda.bat%"
    if exist "!CONDA_BAT!" goto :ACTIVATE
)

:: 2. Check standard installation directories for Anaconda and Miniconda
set "PATHS[0]=%USERPROFILE%\anaconda3"
set "PATHS[1]=%USERPROFILE%\miniconda3"
set "PATHS[2]=%ProgramData%\anaconda3"
set "PATHS[3]=%ProgramData%\miniconda3"
set "PATHS[4]=C:\anaconda3"
set "PATHS[5]=C:\miniconda3"

for /L %%i in (0,1,5) do (
    if exist "!PATHS[%%i]!\condabin\conda.bat" (
        set "CONDA_BAT=!PATHS[%%i]!\condabin\conda.bat"
        goto :ACTIVATE
    )
)

:: 3. Fallback: Check if conda is already in system PATH
where conda >nul 2>&1
if %errorlevel%==0 (
    call conda activate StrainVis_1_4
    goto :RUN
)

echo ERROR: Could not find Conda installation.
echo Please install Anaconda/Miniconda or run this script from an Anaconda Prompt.
pause
exit /b 1

:ACTIVATE
echo Found Conda at: "%CONDA_BAT%"
:: Initialize Conda and activate the target environment
call "%CONDA_BAT%" activate StrainVis_1_4

:RUN
echo.
echo Starting StrainVis...
python run_strainvis.py --port 5005 --show

if %errorlevel% neq 0 (
    echo.
    echo Application stopped or encountered an error.
)

pause