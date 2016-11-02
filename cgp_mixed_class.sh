#!/bin/bash
# cgp_mixed_class.sh

for c in 10 20 50 100 200 500 750
do
  for d in 10 20 50 100 200 250
  do

    echo $c $d
    file_name="all_scaled_C_${c}_D_${d}"

    Rscript GIT/cgp_new_prep.R $c $d $file_name

    export file_name
    sbatch GIT/cgp_mixed_class.cmd
    echo "Done sending mlp job"

  done

done

echo "DONE"
