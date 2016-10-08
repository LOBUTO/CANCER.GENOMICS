#!/bin/bash
# CGP_BASELINE_DL_FEATURES TESTING
# Produces model for most variable cell and drug features
# across entire - target cgp data set

#10 20 50 100 200 300 500 750 900
#10 20 50 100 200
LOG_FOLDER="/tigress/zamalloa/CGP_FILES/CGP_BASELINE_RESULTS/log.Baseline"
for c in 750
do
  for d in 10 50 200
  do
    for drug in Dasatinib Camptothecin Mitomycin_C Imatinib Cisplatin Sunitinib Midostaurin 17-AAG \
    Gefitinib Nilotinib Doxorubicin Vorinostat Gemcitabine Cytarabine \
    Axitinib Vinorelbine Methotrexate Shikonin Etoposide Paclitaxel Embelin Cyclopamine \
    PAC-1 Bleomycin Docetaxel Rapamycin ATRA Sorafenib Erlotinib Temsirolimus Parthenolide Lapatinib \
    Pazopanib Vinblastine Bortezomib Pyrimethamine Elesclomol Roscovitine
    do

      log_file="${LOG_FOLDER}_${drug}_C_${c}_D_${d}.txt"

      if [ ! -f $log_file ]
      then

        echo $drug $c $d
        # Build feature tables
        Rscript GIT/cgp_baseline_mlp.R $drug $c $d
        echo "Done building training sets"

        # Execute mlp model
        export drug c d
        sbatch GIT/cgp_baseline_mlp.cmd
        echo "Done sending mlp job"
      fi

    done
  done
done
