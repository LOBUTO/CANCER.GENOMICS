#!/bin/bash
# cgp_multi_pred.sh

target=$1 #nci_60, ccle, tcga_brca, tcga_luad, tcga_stad, tcga_coad, tcga_all, tcga_multi, self, other
class_mlp=$2 #T/F
met_type=$3 #morgan_bits, morgan_counts
drug_n=$4
cell_n=$5
fusion_n=$6
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, or None
bn_external=$9 #T/F
pca=${10} #T/F
rebalance=${11} #T/F
mf_manual=${12}

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

# zero_PHA-665752 zero_PD-0325901
for samples in zero_PHA-665752 zero_PD-0325901
do
  for c in 950 # Number of cell features
  do
    for d in 0 # Number of drug features
    do
      for r in 16 # Morgan radii settings
      do
        for b in 2048 # Morgan bit settings (Not needed for morgan counts choice)
        do

          echo $c $d  $r $b $samples

          file_tag_1="/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_radii_${r}_bit_${b}"

          if [ "$class_mlp" == "T" ]
          then
            file_tag_2="/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_radii_${r}_bit_${b}"
          else
            file_tag_2="/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS/${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_radii_${r}_bit_${b}"
          fi

          cgp_drug="${file_tag_1}_train_drug"
          cgp_cell="${file_tag_1}_train_cell"
          model_file="${file_tag_2}.pkl"

          if [ "$target" == "tcga" ]
          then
            for cancer in tcga_luad tcga_brca tcga_stad tcga_coad
            do
              echo $cancer
              # Prep data to be predicted on
              Rscript GIT/cgp_multi_pred.R $cancer $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $bn_external $pca $r $b

              script_name="GIT/cgp_multi_pred.py"
              file_name="$cancer $cgp_drug $model_file $class_mlp $bn_external $d"

              export script_name file_name
              sbatch GIT/cgp_mixed_class.cmd
              echo "Done sending multiplicative_fusion predictive mlp job"

            done

          else
            # Prep data to be predicted on
            Rscript GIT/cgp_multi_pred.R $target $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $bn_external $pca $r $b

            script_name="GIT/cgp_multi_pred.py"
            file_name="$target $cgp_drug $model_file $class_mlp $bn_external $d"

            export script_name file_name
            sbatch GIT/cgp_mixed_class.cmd
            echo "Done sending multiplicative_fusion predictive mlp job"
          fi
        done
      done
    done
  done
done
echo "DONE"
