#!/bin/bash
#tcga_model_pred.sh
#Bash script to predict and produce plots for tcga drug predictions on the nci60 based model

# python /home/zamalloa/Documents/FOLDER/GIT/tcga_model_pred_pca.py \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062816_tcga_all_feat_table.csv \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/nci60_all_feat_table_scaled_round.3.csv \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/PCA.MODELS/nci60.pca.3.random208.pc1000.pkl \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/PCA.MODELS/ALL_nci60_pca500_class_model_500.500.pkl \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.CANCER.SAMPLES.csv \
# /home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/nci60_train_scaling500.pkl
#
# Rscript /home/zamalloa/Documents/FOLDER/GIT/tcga_prediction_plot.R

python /tigress/zamalloa/GIT/tcga_model_pred_pca.py \
/tigress/zamalloa/TABLES/TCGA.TRAINING/062816_tcga_all_feat_table.csv \
/tigress/zamalloa/TABLES/TCGA.TRAINING/nci60_all_feat_table_scaled_round.3.csv \
/tigress/zamalloa/TABLES/TCGA.TRAINING/PCA.MODELS/nci60.pca.3.random208.pc1000.pkl \
/tigress/zamalloa/RESULTS/TCGA.TRAINING/ALL_nci60_pca500_class_model_500.500.pkl \
/tigress/zamalloa/TABLES/TCGA.TRAINING/062116.CANCER.SAMPLES.csv \
/tigress/zamalloa/TABLES/TCGA.TRAINING/nci60_train_scaling500.pkl

Rscript /tigress/zamalloa/GIT/tcga_prediction_plot.R

echo "Congrats!!!"
