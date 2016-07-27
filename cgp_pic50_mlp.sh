#!/bin/bash
#cgp_pic50_mlp.sh

# Take arguments
layers=$1

# Execute mlp
# for drug in A-770041 AMG-706 Axitinib AZD-0530 BIBW2992 \
#                    BMS-536924 Bosutinib Erlotinib FTI-277 Imatinib Lapatinib \
#                    NVP-TAE684 OSI-906 Pazopanib PD-173074 PF-02341066
# for drug in PHA-665752 PLX4720 SB590885 Sorafenib Sunitinib WH-4-023 A-443654 \
#                 BX-795 GDC0941 JW-7-52-1 Midostaurin Rapamycin Temsirolimus
for drug in 17-AAG AZD6244 Nilotinib PD-0325901 PD-0332991 Paclitaxel        
do

  echo $drug
  export layers drug

  sbatch /tigress/zamalloa/GIT/cgp_pic50_mlp.cmd

  sleep 1
done

# Exit
echo "submitted to sbatch"
