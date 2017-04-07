#!/bin/bash
# cgp_mixed_class.sh

class_mlp=$1 #F/T
met_type=$2
multiplicative_fusion=$3 #T/F
drug_n=$4 #manual_x_..
cell_n=$5 #manual_x_..
fusion_n=$6 #manual_x_..
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, tcga_coad, tcga_luad, tcga_stad or None
pca=$9 #T/F
rebalance=${10} #F/T
gene_target=${11} #Specific to target dataset: ccle, geeleher_cisplatin, geeleher_docetaxel or None
fold=${12} #fold_all, fold_early, fold_none
mf_manual=${13} #How many mf weights in mf input layer (needs to be used if last layer of drug_n and cell_n are not equal)

export PATH=/usr/local/cuda/bin:/usr/local/cuda/lib64:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

for samples in all
do
  for c in 300 600 # Number of cell features (500 900)
  do
    cn=$c
    ch=$(($cn/2))
    for cell_n in "manual_${cn}"
    do
    for d in 150 290 # Number of drug features (100 290)
    do
      dn=$d
      dh=$(($dn/2))
      for drug_n in "manual_${dn}"
      do
        last_c=${cell_n##m*_}
        last_d=${drug_n##m*_}
        last_total=$(($last_d+$last_c))
        last_half=$(($last_total/2))
      for fusion_n in "manual_${last_total}_${last_total}" "manual_${last_total}_${last_half}"
      do
      for r in 16 # Morgan radii settings
      do
        for b in 2048 # Morgan bit settings (Not needed for morgan counts choice)
        do
          echo $c $cell_n $d $drug_n $fusion_n $samples

          # Prep training sets
          if [ "$multiplicative_fusion" == "T" ]
          then

            if [ "$fold" == "fold_early" ] || [ "$fold" == "fold_none" ]
            then
              echo "$fold"
              file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_gene_target_${gene_target}_fold_${fold}_radii_${r}_bit_${b}"

              Rscript GIT/cgp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $pca $rebalance $r $b $gene_target $fold

              script_name="GIT/cgp_multi_mlp.py"

              file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${fold} ${mf_manual}"

              THEANO_FLAGS='mode=FAST_RUN,allow_gc=True,linker=c,device=gpu,lib.cnmem=1,floatX=float32,nvcc.fastmath=True' python ${script_name} ${file_name}

              echo "Done sending multiplicative_fusion mlp job"
            else
              echo "$fold"
              file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_gene_target_${gene_target}_fold_${fold}_radii_${r}_bit_${b}"

              Rscript GIT/cgp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $pca $rebalance $r $b $gene_target $fold

              for split_fold in 1 2 3 4 5
              do
                echo "$split_fold"

                script_name="GIT/cgp_multi_mlp.py"
                echo "$cell_n"

                file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_gene_target_${gene_target}_fold_${split_fold}_radii_${r}_bit_${b}"
                file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${fold} ${mf_manual}"

                export script_name file_name
                sbatch GIT/cgp_mixed_class.cmd
                echo "Done sending multiplicative_fusion mlp job"
              done
            fi

          else
            Rscript GIT/cgp_new_prep.R $c $d $file_name $met_type $class_mlp $samples $genes

            file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n}"

            if [ "$class_mlp" == "T" ]
            then

              script_name="GIT/cgp_mixed_class.py"

              export script_name file_name
              sbatch GIT/cgp_mixed_class.cmd
              echo "Done sending non-mf class mlp job"

            else

              script_name="GIT/cgp_mixed_reg.py"

              export file_name script_name
              sbatch GIT/cgp_mixed_class.cmd
              echo "Done sending non-mf regression mlp job"
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
echo "DONE"
