#!/bin/bash
# main_pipeline.sh
conda create -c rdkit -n my-rdkit-env rdkit
source activate my-rdkit-env

# Get parameters
cell_array=$1
drug_smiles=$2
target=$3
out_file=$4

# Execute python code
python main_pipeline.py $cell_array $drug_smiles $target $out_file
