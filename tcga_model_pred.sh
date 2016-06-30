#!/bin/bash
#tcga_model_pred.sh
#Bash script to predict and produce plots for tcga drug predictions on the nci60 based model

python /home/zamalloa/Documents/FOLDER/GIT/tcga_model_pred.py \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.TCGA.DRUG.FEAT_SCALED.csv \
/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/ALL_nci60_class_model_.pkl \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.CANCER.SAMPLES.csv

Rscript /home/zamalloa/Documents/FOLDER/GIT/tcga_prediction_plot.R

echo "Congrats!!!"
