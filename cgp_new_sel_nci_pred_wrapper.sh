#!/bin/bash
# cgp_new_sel_nci_pred_wrapper.sh

model_start="CGP_FILES/CGP_NEW_RESULTS/new_cgp_sel_model_new_combat_"
model_end=".pkl"
extra="new_combat"

for drug in Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib

do

  model_file="$model_start$drug$model_end"
  ./GIT/cgp_new_sel_nci_pred.sh $drug $model_file $extra

done

echo "Done with wrapper"
