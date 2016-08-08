#!/bin/bash
# cgp_new_sel_nci_pred.sh
# Function to forecast prediction on nci60 given drug using cgp model

drug=$1
model=$2

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

Rscript GIT/cgp_new_sel_nci_feat.R $drug
echo "Done building training sets"

python GIT/cgp_new_sel_nci_pred.py $drug $model
echo "Done predicting"

Rscript GIT/cgp_new_sel_nci_pred.R $drug
echo "Done graphing"
