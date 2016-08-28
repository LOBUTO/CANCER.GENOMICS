#!/bin/bash
# cgp_new_sel_mlp_tigress.sh

extra=$1
usage=$2
modifier=$3

# for drug in Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
# PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
# Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine

#for drug in PLX4720 TAE684 PHA-665752 Sorafenib PD-0325901 PD-0332991 Lapatinib 17-AAG Erlotinib Nilotinib Paclitaxel

# for drug in Erlotinib Lapatinib Vorinostat Elesclomol ATRA Gefitinib Parthenolide Cyclopamine 17-AAG \
# Sunitinib Gemcitabine Doxorubicin Vinorelbine Vinblastine Mitomycin_C Embelin Paclitaxel
for drug in Erlotinib Lapatinib Gefitinib Sunitinib Gemcitabine Doxorubicin
do

  echo $drug
  Rscript GIT/cgp_new_sel_feat.R $drug $usage $modifier
  echo "Done building training sets"

  export drug usage extra modifier
  sbatch GIT/cgp_new_sel_mlp_tigress.cmd
  echo "Done sending mlp job"
done
