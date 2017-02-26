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

# for samples in zero_17-AAG zero_Nilotinib zero_PD-0325901 zero_PD-0332991 zero_PLX4720 zero_Erlotinib \
# zero_Lapatinib zero_PHA-665752 zero_Paclitaxel zero_Sorafenib zero_TAE684 \
# zero_Crizotinib zero_Nutlin-3 zero_Saracatinib zero_selumetinib
# for samples in percent_all_50 percent_all_20 percent_all_10 percent_all_70
for samples in act_rebalance_5 act_rebalance_4 act_rebalance_3 act_rebalance_2 act_rebalance_1 \
act_rebalance_1.5 act_rebalance_2.5 act_rebalance_3.5 act_rebalance_4.5
do
  for c in 950 # Number of cell features
  do
    for d in 2048 # Number of drug features
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
          # model_file="${file_tag_2}.pkl"

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

            model_files="${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_radii_${r}_bit_${b}"
            model_files=$( ls CGP_FILES/REGRESSION_RESULTS/ | grep $model_files | grep -v combined | grep -v log | grep -v 2048.pkl)

            for model in $model_files
            do

              model_file="CGP_FILES/REGRESSION_RESULTS/${model}"
              echo $model_file

              script_name="GIT/cgp_multi_pred.py"
              file_name="$target $cgp_drug $model_file $class_mlp $bn_external $d"

              export script_name file_name
              sbatch GIT/cgp_mixed_class.cmd
              echo "Done sending multiplicative_fusion predictive mlp job"
            done
          fi
        done
      done
    done
  done
done
echo "DONE"
