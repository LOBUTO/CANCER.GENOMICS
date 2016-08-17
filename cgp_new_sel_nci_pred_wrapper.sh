#!/bin/bash
# cgp_new_sel_nci_pred_wrapper.sh

model_start="CGP_FILES/CGP_NEW_RESULTS/new_cgp_sel_model_new_combat_"
model_end=".pkl"
extra="new_combat"
usage=$1

for drug in Dasatinib Camptothecin Mitomycin_C Imatinib Cisplatin Sunitinib Midostaurin 17-AAG \
Gefitinib Nilotinib Doxorubicin Vorinostat Gemcitabine Cytarabine \
Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine
# for drug in PLX4720 TAE684 PHA-665752 Sorafenib PD-0325901 PD-0332991 Lapatinib 17-AAG Erlotinib Nilotinib Paclitaxel
do

  model_file="$model_start$drug$model_end"
  ./GIT/cgp_new_sel_nci_pred.sh $drug $model_file $extra $usage

done

echo "Done with wrapper"
