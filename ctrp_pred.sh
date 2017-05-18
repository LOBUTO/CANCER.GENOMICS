#!/bin/bash
# ctrp_pred.sh

target=$1 #nci_60, ccle, tcga_brca, tcga_luad, tcga_stad, tcga_coad, tcga_all, tcga_multi, self, other
class_mlp=$2 #T/F
met_type=$3 #morgan_bits, morgan_counts
drug_n=$4
cell_n=$5
fusion_n=$6
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, or None
genespca=$9 #T/F
drugspca=${10}
fold=${11}
mf_manual=${12}

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

for samples in all
do
  for c in 800 # Number of cell features
  do
    cn=$c
    ch=$(($cn/2))
    cnh=$(($cn+$ch))
    for cell_n in manual_800_400_40 manual_800_400_20 manual_800_200_40 manual_800_200_20 #"manual_${cnh}" # "manual_${cn}_${ch}"
    do
    for d in 512 # Number of drug features
    do
      dn=$d
      dh=$(($dn/2))
      dnh=$(($dn+$dh))
      for drug_n in manual_512_200_40 manual_512_200_20  #"manual_${dnh}" # "manual_${dn}_${dh}"
      do
        last_c=${cell_n##m*_}
        last_d=${drug_n##m*_}
        last_total=$(($last_d+$last_c))
        last_half=$(($last_total/2))
        last_total_half=$(($last_total+$last_half))
        last_multi=$(($last_d*$last_c))
        last_multi_half=$(($last_multi/2))

      for fusion_n in "manual_${last_multi}" "manual_${last_multi}_${last_multi_half}" "manual_${last_multi_half}_${last_multi_half}" #"manual_${last_total}_${last_half}" "manual_${last_total}_${last_half}_${last_half}" "manual_${last_total_half}_${last_total}" "manual_${last_total}" "manual_${last_total}_${last_total}" "manual_${last_total}_${last_total}_${last_total}"
      do
      for r in 2 # Morgan radii settings
      do
        for b in 512 # Morgan bit settings (Not needed for morgan counts choice)
        do

          # for th_split in th_1.1_0.4 th_0.9_0.5 th_0.5_0.4 th_1.1_0.3 th_1.1_0.2 th_1.1_0.1
          # for th_split in th_1.1_0.4 th_1.1_0.3 th_1.1_0.2 th_1.1_0.1 th_1.0_0.3 th_1.0_0.2 th_1.0_0.1 th_0.8_0.4 th_0.8_0.5 \
          # th_0.8_0.6 th_1.2_0.6 th_1.2_0.7 th_1.2_0.8
          for th_split in th_0_0
          do
            echo $c $cell_n $d $drug_n $fusion_n $samples $th_split

            file_tag_1="/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"

            if [ "$class_mlp" == "T" ]
            then
              results_folder="CLASS_RESULTS"
            else
              results_folder="REGRESSION_RESULTS"
            fi
            file_tag_2="/tigress/zamalloa/CGP_FILES/${results_folder}/${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"

            cgp_drug="${file_tag_1}_train_drug"
            cgp_cell="${file_tag_1}_train_cell"

            if [ "$target" == "tcga" ]
            then
              for cancer in tcga_luad tcga_brca tcga_stad tcga_coad
              do
                echo $cancer
                # Prep data to be predicted on
                Rscript GIT/ctrp_pred.R $cancer $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $genespca $drugspca $r $b $genes

                model_files="${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"
                model_files=$( ls CGP_FILES/${results_folder}/ | grep $model_files | grep -v combined | grep -v log | grep -v 512.pkl)

                for model in $model_files
                do

                  exists=${model/.pkl/}
                  exists="${exists}_${cancer}_PREDICTION"
                  exists="/tigress/zamalloa/PREDICTIONS/${exists}"
                  echo $exists

                  if [ ! -f "$exists" ]
                  then
                    model_file="CGP_FILES/${results_folder}/${model}"
                    echo $model_file

                    script_name="GIT/cgp_multi_pred.py"
                    file_name="$cancer $cgp_drug $model_file $class_mlp $d"

                    export script_name file_name
                    sbatch GIT/cgp_mixed_class.cmd
                    echo "Done sending multiplicative_fusion predictive mlp job"
                  fi
                done
              done

            else
              # Prep data to be predicted on
              Rscript GIT/ctrp_pred.R $target $met_type $cgp_drug $cgp_cell $class_mlp $batch_norm $genespca $drugspca $r $b $genes

              model_files="${samples}_scaled_C_${c}_${mm}_${d}_mf_T_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"
              model_files=$( ls CGP_FILES/${results_folder}/ | grep $model_files | grep -v combined | grep -v log | grep -v 512.pkl)

              for model in $model_files
              do

                exists=${model/.pkl/}
                exists="${exists}_${target}_PREDICTION"
                exists="/tigress/zamalloa/PREDICTIONS/${exists}"
                echo $exists

                if [ ! -f "$exists" ]
                then
                  model_file="CGP_FILES/${results_folder}/${model}"
                  echo $model_file

                  script_name="GIT/cgp_multi_pred.py"
                  file_name="$target $cgp_drug $model_file $class_mlp $d"

                  export script_name file_name
                  sbatch GIT/cgp_mixed_class.cmd
                  echo "Done sending multiplicative_fusion predictive mlp job"
                fi
              done
            fi
          done
        done
      done
      done
      done
    done
    done
  done
done
echo "DONE"
