#!/bin/bash
#cgp.cor.AUC_class_pred.sh
#Bash to predict and produce plot for cgp.cor.AUC based on the nci60 based model

python /home/zamalloa/Documents/FOLDER/GIT/cgp_model_pred_pca.py \
/home/zamalloa/Documents/FOLDER/CGP_FILES/cgp_all_feat_class \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/nci60_all_feat_table_scaled_round.3.csv \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/PCA.MODELS/nci60.pca.3.random208.pc1000.pkl \
/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/PCA.MODELS/ALL_nci60_pca700_class_model_700.700.pkl

Rscript /home/zamalloa/Documents/FOLDER/GIT/cgp.cor.AUC_class_pred_plot.R

echo "Congrats!!"
