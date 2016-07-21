#!/bin/bash
#cgp_pic50_mlp.sh

# Take arguments
drug=$1
layers=$2

# Execute mlp
export layers drug

sbatch /tigress/zamalloa/GIT/cgp_pic50_mlp.cmd

# Exit
echo "submitted to sbatch"
