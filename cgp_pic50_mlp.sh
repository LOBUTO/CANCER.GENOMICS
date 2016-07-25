#!/bin/bash
#cgp_pic50_mlp.sh

# Take arguments
layers=$1

# Execute mlp
for drug in A-443654 A-770041 AMG-706 Axitinib AZD-0530 BIBW2992 \
                   BMS-536924 Bosutinib Erlotinib FTI-277 Imatinib Lapatinib \
                   NVP-TAE684 OSI-906 Pazopanib PD-173074 PF-02341066
do

  echo $drug
  export layers drug

  sbatch /tigress/zamalloa/GIT/cgp_pic50_mlp.cmd

  sleep 1
done

# Exit
echo "submitted to sbatch"
