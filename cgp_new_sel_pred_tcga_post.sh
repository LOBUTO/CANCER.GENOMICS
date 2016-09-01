#!/bin/bash
# cgp_new_sel_pred_tcga_post.sh

modifier=$1
extra=$2

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

for drug in Anastrozole Docetaxel Paclitaxel Tamoxifen
do

    cancer="brca"
    echo $cancer

    python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

    Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

for drug in Paclitaxel Carboplatin
do

  cancer="ucec"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

for drug in Alimta Paclitaxel Carboplatin Cisplatin
do

  cancer="luad"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

echo "DONE"
