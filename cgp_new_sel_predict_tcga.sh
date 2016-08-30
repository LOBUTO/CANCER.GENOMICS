#!/bin/bash
# cgp_new_sel_predict_tcga.sh

cancer=$1
cancer_exp=$2
target_drug=$3

usage="tcga_$cancer"
modifier=$4

# Set up entry table
Rscript GIT/cgp_new_sel_predict_tcga_prep.R $cancer $cancer_exp $target_drug
echo "Done building cgp feat table"

# Split for training
Rscript GIT/cgp_new_sel_feat.R $target_drug $usage $modifier
