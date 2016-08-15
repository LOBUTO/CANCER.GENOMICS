#!/bin/bash
# cgp_new_sel_mlp_tigress.sh

drug=$1
extra=$2

for drug in Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine

do
  Rscript GIT/cgp_new_sel_feat.R $drug
  echo "Done building training sets"

  export drug extra
  sbatch GIT/cgp_new_sel_mlp_tigress.cmd $drug $extra
  echo "Done sending mlp job"
done
