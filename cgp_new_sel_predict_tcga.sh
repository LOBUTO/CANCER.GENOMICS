#!/bin/bash
# cgp_new_sel_predict_tcga.sh

function tcga {

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

}

modifier=$1
extra=$2

# cancer="COAD"
# for i in CAPECITABINE
# do
#   tcga $cancer $i $modifier $extra
# done

cancer="COADREAD"
for i in CAPECITABINE
do
  tcga $cancer $i $modifier $extra
done

cancer="LUAD"
for i in Paclitaxel
do
  tcga $cancer $i $modifier $extra
done

cancer="LUSC"
for i in Carboplatin Gemcitabine
do
  tcga $cancer $i $modifier $extra
done

cancer="MESO"
for i in Cisplatin
do
  tcga $cancer $i $modifier $extra
done

cancer="PAAD"
for i in Gemcitabine 5-Fluorouracil
do
  tcga $cancer $i $modifier $extra
done

cancer="READ"
for i in Leucovorin 5-Fluorouracil
do
  tcga $cancer $i $modifier $extra
done

cancer="READ"
for i in Leucovorin 5-Fluorouracil
do
  tcga $cancer $i $modifier $extra
done

cancer="STAD"
for i in Leucovorin 5-Fluorouracil Cisplatin Epirubicin
do
  tcga $cancer $i $modifier $extra
done

cancer="STES"
for i in 5-Fluorouracil CAPECITABINE Cisplatin Epirubicin Etoposide Leucovorin
do
  tcga $cancer $i $modifier $extra
done

cancer="UCEC"
for i in Carboplatin Paclitaxel
do
  tcga $cancer $i $modifier $extra
done

cancer="UCS"
for i in Carboplatin
do
  tcga $cancer $i $modifier $extra
done
