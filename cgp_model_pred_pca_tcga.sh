#!/bin/bash
#cgp.cor.AUC_class_pred.sh
#Bash to predict and produce plot for cgp.cor.AUC based on the nci60 based model

for p in 500 800 1000
do
  python /home/zamalloa/Documents/FOLDER/GIT/cgp_model_pred_pca_tcga.py \
  /home/zamalloa/Documents/FOLDER/CGP_FILES/cgp_all_feat \
  /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062816_tcga_all_feat_table_rescaled.csv \
  /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/PCA.MODELS/cgp_pca_1000.pkl \
  /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.CANCER.SAMPLES.csv \
  $p
done

Rscript /home/zamalloa/Documents/FOLDER/GIT/cgp_model_tcga_plot.R

echo "Congrats!!"
