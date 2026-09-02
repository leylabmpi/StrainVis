@echo off
echo Creating Conda environment. This may take a few minutes...
call conda env create -f strainvis.yml
echo Setup complete! You can close this window.
pause