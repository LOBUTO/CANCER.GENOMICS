#!/bin/bash
# ccle_pca_predict.sh

drug=$1
m_pca=$2
g_pca=$3

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

Rscript GIT/ccle_pca_predict_preprocess.R $drug $m_pca $g_pca
echo "Done preprocessing"

python GIT/ccle_pca_predict.py $drug $m_pca $g_pca
echo "Done building predictions"

Rscript GIT/ccle_pca_predict.R $drug $m_pca $g_pca
echo "Done plotting predictions"
