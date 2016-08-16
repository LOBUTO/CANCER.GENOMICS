#!/bin/bash
# cgp_new_sel_mlp_tigress.sh

extra=$1
usage=$2

# for drug in Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
# PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
# Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine
for drug in PLX4720 TAE684 PHA-665752 Sorafenib PD-0325901 PD-0332991 Lapatinib 17-AAG Erlotinib Nilotinib Paclitaxel
do
  Rscript GIT/cgp_new_sel_feat.R $drug $usage
  echo "Done building training sets"

  export drug usage extra
  sbatch GIT/cgp_new_sel_mlp_tigress.cmd
  echo "Done sending mlp job"
done
