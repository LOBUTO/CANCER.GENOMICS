#!/bin/bash
#tcga_model_pred.sh
#Bash script to predict and produce plots for tcga drug predictions on the nci60 based model

python /home/zamalloa/Documents/FOLDER/GIT/tcga_model_pred_pca.py \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062816_tcga_all_feat_table.csv \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/nci60_all_feat_table_scaled_round.2.csv \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/PCA.MODELS/nci60.pca.random214.pc500.pkl \
/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/ALL_nci60_pca_class_model_300.300.pkl \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.CANCER.SAMPLES.csv

Rscript /home/zamalloa/Documents/FOLDER/GIT/tcga_prediction_plot.R

echo "Congrats!!!"
