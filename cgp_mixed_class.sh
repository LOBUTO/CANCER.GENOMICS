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

# for samples in FK866 IPA-3 NSC-207895 UNC0638 CX-5461 Trametinib SNX-2112 OSI-027 \
# QS11 AT-7519 PAC-1 SN-38 PI-103 I-BET-762 5-Fluorouracil PHA-793887 YM201636 \
# LY317615 TAK-715 RDEA119 Gemcitabine Bleomycin CHIR-99021 VX-11e EKB-569 GSK-650394 \
# 17-AAG AG-014699 GDC0941 Y-39983
# for samples in Gemcitabine Sunitinib Doxorubicin Vinorelbine Vinblastine Midostaurin Paclitaxel \
# Camptothecin Embelin Bleomycin Axitinib Docetaxel Nilotinib Sorafenib Cytarabine \
# Shikonin Roscovitine Etoposide Pyrimethamine Methotrexate PAC-1 Temsirolimus Rapamycin Bortezomib \
# Imatinib Pazopanib Dasatinib Lapatinib Erlotinib Vorinostat Cyclopamine Cisplatin \
# Elesclomol 17-AAG ATRA Gefitinib Parthenolide

if [ "$met_type" == "morgan_bits" ]
then
  mm="MB"
else
  mm="MC"
fi

# for samples in zero_17-AAG zero_Nilotinib zero_PD-0325901 zero_PD-0332991 zero_PLX4720 zero_Erlotinib \
# zero_Lapatinib zero_PHA-665752 zero_Paclitaxel zero_Sorafenib zero_TAE684 \
# zero_Crizotinib zero_Nutlin-3 zero_Saracatinib zero_selumetinib
# for samples in zero_5-Fluorouracil zero_Cisplatin zero_Gemcitabine zero_Temozolomide
#for samples in percent_all_50 percent_all_20 percent_all_10 percent_all_70
# for samples in act_rebalance_5 act_rebalance_4 act_rebalance_3 act_rebalance_2 act_rebalance_1 \
# act_rebalance_1.5 act_rebalance_2.5 act_rebalance_3.5 act_rebalance_4.5
# for samples in act_rebalancetop_10 act_rebalancetop_20 act_rebalancetop_40 act_rebalancetop_50 \
# act_rebalancetop_60 act_rebalancetop_80 act_rebalancetop_100
# for samples in zero_Docetaxel zero_Cisplatin zero_Bortezomib zero_Erlotinib
# for cell_n in manual_950 manual_950_950 manual_50_50 manual_50_50_50 manual_50 manual_100 manual_100_100 manual_100_100_100 \
# manual_200 manual_200_200 manual_200_200_200 manual_500 manual_500_500 manual_500_500_500 manual_750 manual_750_750_750
# for cell_n in manual_500 manual_500_500 manual_50_50 manual_50_50_50 manual_50 manual_100 manual_100_100 manual_100_100_100 \
# manual_200 manual_200_200 manual_200_200_200
for samples in all
do
  for c in 300 600 950 # Number of cell features mb_early_none(300 600 950) folds(100 300 600 950)
  do
    cn=$c
    ch=$(($cn/2))
    cnh=$(($cn+$ch))
    for cell_n in "manual_${cnh}" # "manual_${cn}_${ch}"
    do
    for d in 50 100 290 # Number of drug features mb_early_none(50 100 290) folds(50 150 290)
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
          echo $c $cell_n $d $drug_n $fusion_n $samples

          # Prep training sets
          if [ "$multiplicative_fusion" == "T" ]
          then

            if [ "$fold" == "fold_early" ] || [ "$fold" == "fold_none" ]
            then
              echo "$fold"
              file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_rebalance_${rebalance}_gene_target_${gene_target}_fold_${fold}_radii_${r}_bit_${b}"

              Rscript GIT/cgp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $pca $rebalance $r $b $gene_target $fold

              script_name="GIT/cgp_multi_mlp_2.py"

              file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${fold} ${mf_manual}"

              export script_name file_name
              sbatch GIT/cgp_mixed_class.cmd
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
