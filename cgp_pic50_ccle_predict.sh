#!/bin/bash
# cgp_pic50_ccle_predict.sh
# Master wrapper for python script to predict ccle drug activity using cgp models
# DEPRECATED as python script was self sufficient

module load python
module load cudatoolkit

for drug in 17-AAG AEW541 AZD0530 AZD6244 Erlotinib Irinotecan L-685458 \
              LBW242 Lapatinib Nilotinib Nutlin-3 PD-0325901 PD-0332991 \
              PF2341066 PHA-665752 PLX4720 Paclitaxel Panobinostat RAF265 \
              Sorafenib TAE684 TKI258 Topotecan ZD-6474
do

  echo $drug
  python /tigress/zamalloa/GIT/cgp_pic50_ccle_predict.py $drug

done
