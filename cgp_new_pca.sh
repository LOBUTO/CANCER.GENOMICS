#!/bin/bash
# cgp_new_pca.sh

drug=$1
m_pca=$2
g_pca=$3

cgp_new_pca_feat.R $drug $m_pca $g_pca
echo "Done building training sets"

cgp_new_pca_mlp.py $drug $m_pca $g_pca
echo "Done sending mlp job"
