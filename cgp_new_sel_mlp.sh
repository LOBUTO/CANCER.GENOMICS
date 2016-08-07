#!/bin/bash
# cgp_new_sel_mlp.sh

drug=$1

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

Rscript GIT/cgp_new_sel_feat.R $drug
echo "Done building training sets"

python GIT/cgp_new_sel_mlp.py $drug
echo "Done sending mlp job"
