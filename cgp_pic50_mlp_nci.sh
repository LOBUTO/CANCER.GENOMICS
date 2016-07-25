#!/bin/bash
#cgp_pic50_mlp_nci.sh

# Take arguments
layers=$1

# Execute mlp
for drug in 718781
do

  echo $drug
  export layers drug

  sbatch /tigress/zamalloa/GIT/cgp_pic50_mlp.cmd

  sleep 1
done

# Exit
echo "submitted to sbatch"
