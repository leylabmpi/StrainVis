#!/bin/bash
cd "$(dirname "$0")"

echo "Starting StrainVis..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate StrainVis_1_4

python run_strainvis.py --port 5005 --show