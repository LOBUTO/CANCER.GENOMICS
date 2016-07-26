#!/bin/bash
#cgp_pic50_mlp_ccle.sh

# Take arguments
layers=$1

# Execute mlp
for drug in 17-AAG AEW541 AZD0530 AZD6244 Erlotinib Irinotecan L-685458 \
              LBW242 Lapatinib Nilotinib Nutlin-3 PD-0325901 PD-0332991 \
              PF2341066 PHA-665752 PLX4720 Paclitaxel Panobinostat RAF265 \
              Sorafenib TAE684 TKI258 Topotecan ZD-6474
do

  echo $drug
  export layers drug

  sbatch /tigress/zamalloa/GIT/cgp_pic50_mlp_ccle.cmd

  sleep 1
done

# Exit
echo "submitted to sbatch"
