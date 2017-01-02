#!/bin/bash
# cgp_multi_pred.sh

target=$1 #nci_60, ccle, tcga_brca, tcga_luad, tcga_stad, tcga_coad, tcga(for all cancers), self, other
class_mlp=$2 #T/F
met_type=$3 #morgan_bits, morgan_counts
drug_n=$4
cell_n=$5
fusion_n=$6
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, or None
bn_external=$9 #T/F
mf_manual=${10}

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

for samples in semi_split
do
  for cell in 500 950
  do
    for drug in 500 1000
    do

      echo $cell $drug $samples

      file_tag_1="/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/${samples}_scaled_C_${cell}_${mm}_${drug}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}"
      file_tag_2="/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/${samples}_scaled_C_${cell}_${mm}_${drug}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}"

      cgp_drug="${file_tag_1}_train_drug"
      cgp_cell="${file_tag_1}_train_cell"
      model_file="${file_tag_2}.pkl"

      if [ "$target" == "tcga" ]
      then
        for cancer in tcga_luad tcga_brca tcga_stad tcga_coad
        do
          echo $cancer
          # Prep data to be predicted on
          Rscript GIT/cgp_multi_pred.R $cancer $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $bn_external

          module load cudatoolkit
          module load python

          THEANO_FLAGS='mode=FAST_RUN,allow_gc=False,linker=c,device=gpu,lib.cnmem=1,floatX=float32,nvcc.fastmath=True' python GIT/cgp_multi_pred.py $cancer $cgp_drug $model_file $class_mlp $bn_external
        done

      else
        # Prep data to be predicted on
        Rscript GIT/cgp_multi_pred.R $target $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $bn_external

        module load cudatoolkit
        module load python

        THEANO_FLAGS='mode=FAST_RUN,allow_gc=False,linker=c,device=gpu,lib.cnmem=1,floatX=float32,nvcc.fastmath=True' python GIT/cgp_multi_pred.py $target $cgp_drug $model_file $class_mlp $bn_external
      fi

    done
  done
done
echo "DONE"
