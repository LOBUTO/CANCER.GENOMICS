#!/bin/bash
# cgp_new_sel_predict_tcga.sh

cancer=$1
drug=$2

usage="tcga_$cancer"
modifier=$3
extra=$4

# Set up entry table
Rscript GIT/cgp_new_sel_predict_tcga_prep.R $cancer $drug
echo "Done building cgp feat table"

# Split for training
feat_table="/tigress/zamalloa/CGP_FILES/${cancer}_all_cgp_new_${drug}.rds"
Rscript GIT/cgp_new_sel_feat.R $drug $usage $modifier $feat_table
echo "Done splitting drug table for training"

# Train
export drug usage extra modifier
sbatch GIT/cgp_new_sel_mlp_tigress.cmd
echo "Done sending mlp job"
