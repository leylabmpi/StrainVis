@echo off
echo Starting StrainVis...
call conda activate StrainVis_1_4
python run_strainvis.py --port 5005 --show
pause