#!/bin/bash
#cgp_pic50_train.sh

#for drug in 17-AAG AZD6244 Nilotinib PD-0325901 PD-0332991 PLX4720 Erlotinib Lapatinib PHA-665752 Paclitaxel Sorafenib
for drug in 17-AAG AZD6244 Nilotinib PD-0325901 PD-0332991 Paclitaxel
do

  echo $drug
  Rscript /tigress/zamalloa/GIT/cgp_pic50_train.R $drug

  module load cudatoolkit
  module load python
  python /tigress/zamalloa/GIT/cgp_pic50_train.py $drug

done

echo "Done"
