#!/bin/bash
# cgp_mixed_class.sh

met_type=="morgan_bits"
for c in 10 20 50 100 200 500 750
do
  for d in 32 64 128 256 512 1024 2048
  do

    echo $c $d
    file_name="all_scaled_C_${c}_MB_${d}"

    Rscript GIT/cgp_new_prep.R $c $d $file_name $met_type

    export file_name
    sbatch GIT/cgp_mixed_class.cmd
    echo "Done sending mlp job"

  done

done

echo "DONE"
