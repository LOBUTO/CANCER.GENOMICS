#!/bin/bash
# cgp_new_sel_pred_tcga_post.sh

cancer=$1
modifier=$2
extra=$3

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

for drug in Docetaxel Paclitaxel Tamoxifen
do

  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done
echo "DONE"
