#!/bin/bash
#cgp.cor.AUC_class_pred.sh
#Bash to predict and produce plot for cgp.cor.AUC based on the nci60 based model

python /home/zamalloa/Documents/FOLDER/GIT/cgp.cor.AUC_class_pred.py \
/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/061716.CGP.COR.AUC.CLASS.csv \
/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/ALL_nci60_class_model_.pkl

Rscript /home/zamalloa/Documents/FOLDER/GIT/cgp.cor.AUC_class_pred_plot.R

echo "Congrats!!"
