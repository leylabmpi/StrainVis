#!/bin/bash
cd "$(dirname "$0")"

echo "Creating Conda environment. This may take a few minutes..."
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env create -f strainvis.yml
echo "Setup complete! You can close this window."