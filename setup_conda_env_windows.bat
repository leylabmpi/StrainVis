@echo off
setlocal enabledelayedexpansion

:: Change directory to the folder where the batch script resides
cd /d "%~dp0"

echo Searching for Conda installation...

:: 1. Check if conda is already in PATH or CONDA_EXE is set
if defined CONDA_EXE (
    set "CONDA_BAT=%CONDA_EXE:\Scripts\conda.exe=\condabin\conda.bat%"
    if exist "!CONDA_BAT!" goto :FOUND
)

:: 2. Check standard installation paths
set "PATHS[0]=%USERPROFILE%\anaconda3"
set "PATHS[1]=%USERPROFILE%\miniconda3"
set "PATHS[2]=%ProgramData%\anaconda3"
set "PATHS[3]=%ProgramData%\miniconda3"
set "PATHS[4]=C:\anaconda3"
set "PATHS[5]=C:\miniconda3"

for /L %%i in (0,1,5) do (
    if exist "!PATHS[%%i]!\condabin\conda.bat" (
        set "CONDA_BAT=!PATHS[%%i]!\condabin\conda.bat"
        goto :FOUND
    )
)

echo ERROR: Could not find Conda installation.
echo Please run this script from an Anaconda Prompt.
pause
exit /b 1

:FOUND
echo Found Conda at: "%CONDA_BAT%"
echo Creating Conda environment. This may take a few minutes...

:: Activate Conda hook and create environment
call "%CONDA_BAT%" env create -f strainvis.yml

echo Setup complete! You can close this window.
pause