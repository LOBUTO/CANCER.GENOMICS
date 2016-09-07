#!/bin/bash
# cgp_new_sel_pred_tcga_post.sh

modifier=$1
extra=$2

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

for drug in Carboplatin Cisplatin Gemcitabine
do

    cancer="BLCA"
    echo $cancer

    python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

    Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

for drug in Cyclophosphamide Paclitaxel Tamoxifen
do

    cancer="BRCA"
    echo $cancer

    python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

    Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

for drug in Cisplatin
do

  cancer="CESC"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in 5-Fluorouracil CAPECITABINE Leucovorin MFCD00866327
do

  cancer="COAD"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

echo "DONE"


for drug in 5-Fluorouracil CAPECITABINE Leucovorin Irinotecan MFCD00866327
do

  cancer="COADREAD"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Temozolomide
do

  cancer="GBMLGG"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Temozolomide
do

  cancer="LGG"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Alimta Carboplatin Cisplatin Paclitaxel
do

  cancer="LUAD"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Carboplatin Gemcitabine
do

  cancer="LUSC"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Cisplatin
do

  cancer="MESO"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in 5-Fluorouracil Gemcitabine
do

  cancer="PAAD"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in 5-Fluorouracil Leucovorin
do

  cancer="READ"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in 5-Fluorouracil Cisplatin Epirubicin Leucovorin
do

  cancer="STAD"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in 5-Fluorouracil CAPECITABINE Cisplatin Epirubicin Etoposide Leucovorin
do

  cancer="STES"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Carboplatin Paclitaxel
do

  cancer="UCEC"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done


for drug in Carboplatin
do

  cancer="UCS"
  echo $drug
  python GIT/cgp_new_sel_pred_tcga_post.py $drug $cancer $modifier $extra

  Rscript GIT/cgp_new_sel_pred_tcga_plot.R $drug $cancer $extra

done

echo "DONE"
