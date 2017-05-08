#!/bin/bash
# ctrp.sh

class_mlp=$1 #F/T
met_type=$2
multiplicative_fusion=$3 #T/F
drug_n=$4 #manual_x_..
cell_n=$5 #manual_x_..
fusion_n=$6 #manual_x_..
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, tcga_coad, tcga_luad, tcga_stad or None
genespca=$9 #T/F
drugspca=${10} #T/F
fold=${11} #fold_all, fold_early, fold_none
mf_manual=${12} #How many mf weights in mf input layer (needs to be used if last layer of drug_n and cell_n are not equal)

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

for samples in all
do
  for c in 400 800
  do
    cn=$c
    ch=$(($cn/2))
    cnh=$(($cn+$ch))
    for cell_n in "manual_${cnh}" # "manual_${cn}_${ch}"
    do
    for d in 480
    do
      dn=$d
      dh=$(($dn/2))
      dnh=$(($dn+$dh))
      for drug_n in "manual_${dnh}" # "manual_${dn}_${dh}"
      do
        last_c=${cell_n##m*_}
        last_d=${drug_n##m*_}
        last_total=$(($last_d+$last_c))
        last_half=$(($last_total/2))
        last_total_half=$(($last_total+$last_half))
      for fusion_n in "manual_${last_total}_${last_half}" # "manual_${last_total}_${last_half}_${last_half}" "manual_${last_total_half}_${last_total}" # "manual_${last_total}" "manual_${last_total}_${last_total}" "manual_${last_total}_${last_total}_${last_total}"
      do
      for r in 16 # Morgan radii settings - 16
      do
        for b in 2048 # Morgan bit settings (Not needed for morgan counts choice) - 2048
        do

          for th_split in th_1.1_0.4 th_0.9_0.5 th_0.5_0.4 th_1.1_0.3 th_1.1_0.2 th_1.1_0.1
          do
          echo $c $cell_n $d $drug_n $fusion_n $samples

          # Prep training sets
          if [ "$multiplicative_fusion" == "T" ]
          then

            if [ "$fold" == "fold_early" ] || [ "$fold" == "fold_none" ]
            then
              echo "$fold"
              file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"

              Rscript GIT/ctrp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $genespca $drugspca $r $b $fold $th_split

              script_name="GIT/ctrp_multi_mlp.py"

              file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${fold} ${mf_manual}"

              export script_name file_name
              sbatch GIT/cgp_mixed_class.cmd
              echo "Done sending multiplicative_fusion mlp job"
            else
              echo "$fold"
              file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${fold}_thsplit_${th_split}_radii_${r}_bit_${b}"

              Rscript GIT/ctrp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $genespca $drugspca $r $b $fold $th_split

              for split_fold in 1 2 3 4 5
              do
                echo "$split_fold"

                script_name="GIT/ctrp_multi_mlp.py"
                echo "$cell_n"

                file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_genespca_${genespca}_drugspca_${drugspca}_fold_${split_fold}_thsplit_${th_split}_radii_${r}_bit_${b}"
                file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${fold} ${mf_manual}"

                export script_name file_name
                sbatch GIT/cgp_mixed_class.cmd
                echo "Done sending multiplicative_fusion mlp job"
              done
            fi
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
