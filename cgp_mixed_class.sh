#!/bin/bash
# cgp_mixed_class.sh

class_mlp=$1
met_type=$2
multiplicative_fusion=$3
drug_n=$4
cell_n=$5
fusion_n=$6 #How many fusion neurons in hidden layer
genes=$7 #F/T
batch_norm=$8 #cgp_nci60, cgp_ccle, tcga_brca, tcga_coad, tcga_luad, tcga_stad or None
pca=$9
mf_manual=${10} #How many mf weights in mf input layer (needs to be used if last layer of drug_n and cell_n are not equal)

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

for samples in semi_split
do
  for c in 10 25 50 100 # Number of cell features
  do
    for d in 10 20 50 100 # Number of drug features
    do
      for r in 1 2 4 8 12 16 # Morgan radii settings
      do
        for b in 256 512 1024 2048 # Morgan bit settings (Not needed for morgan counts choice)
        do

          echo $c $d $samples
          file_name="${samples}_scaled_C_${c}_${mm}_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}_bn_${batch_norm}_pca_${pca}_radii_${r}_bit_${b}"

          # Prep training sets
          if [ "$multiplicative_fusion" == "T" ]
          then
            Rscript GIT/cgp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes $batch_norm $pca $r $b

            script_name="GIT/cgp_multi_mlp.py"

            file_name="${file_name} ${drug_n} ${cell_n} ${fusion_n} ${class_mlp} ${d} ${c} ${mf_manual}"

            export script_name file_name
            sbatch GIT/cgp_mixed_class.cmd
            echo "Done sending multiplicative_fusion mlp job"

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
echo "DONE"
