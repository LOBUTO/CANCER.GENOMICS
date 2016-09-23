#!/bin/bash
# CGP_BASELINE_DL_FEATURES TESTING
# Produces model for most variable cell and drug features
# across entire - target cgp data set

for c in 10 20 50 100 200
do
  for d in 10 20 50 100 200
  do
    for drug in Dasatinib Camptothecin Mitomycin_C Imatinib Cisplatin Sunitinib Midostaurin 17-AAG \
    Gefitinib Nilotinib Doxorubicin Vorinostat Gemcitabine Cytarabine \
    Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
    PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
    Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine
    do

      echo $drug $c $d
      # Build feature tables
      cgp_baseline_mlp.R $drug $c $d
      echo "Done building training sets"

      # Execute mlp model
      export drug c d
      cgp_baseline_mlp.cmd
      echo "Done sending mlp job"

    done
  done
done
