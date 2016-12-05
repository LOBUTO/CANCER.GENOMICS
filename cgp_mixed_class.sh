#!/bin/bash
# cgp_mixed_class.sh

class_mlp=$1
met_type=$2
multiplicative_fusion=$3
drug_n=$4
cell_n=$5
fusion_n=$6
genes=$7
mf_manual=$8

# for samples in FK866 IPA-3 NSC-207895 UNC0638 CX-5461 Trametinib SNX-2112 OSI-027 \
# QS11 AT-7519 PAC-1 SN-38 PI-103 I-BET-762 5-Fluorouracil PHA-793887 YM201636 \
# LY317615 TAK-715 RDEA119 Gemcitabine Bleomycin CHIR-99021 VX-11e EKB-569 GSK-650394 \
# 17-AAG AG-014699 GDC0941 Y-39983
# for samples in Gemcitabine Sunitinib Doxorubicin Vinorelbine Vinblastine Midostaurin Paclitaxel \
# Camptothecin Embelin Bleomycin Axitinib Docetaxel Nilotinib Sorafenib Cytarabine \
# Shikonin Roscovitine Etoposide Pyrimethamine Methotrexate PAC-1 Temsirolimus Rapamycin Bortezomib \
# Imatinib Pazopanib Dasatinib Lapatinib Erlotinib Vorinostat Cyclopamine Cisplatin \
# Elesclomol 17-AAG ATRA Gefitinib Parthenolide
for samples in all
do
  for c in 200
  do
    for d in 0
    do

      echo $c $d $samples
      file_name="${samples}_scaled_C_${c}_MB_${d}_mf_${multiplicative_fusion}_dn_${drug_n}_cn_${cell_n}_fn_${fusion_n}_mf_manual_${mf_manual}_genes_${genes}"

      # Prep training sets
      if [ "$multiplicative_fusion" == "T" ]
      then
        Rscript GIT/cgp_new_prep_mf.R $c $d $file_name $met_type $class_mlp $samples $genes

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
echo "DONE"
