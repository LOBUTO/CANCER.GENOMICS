#!/bin/bash
#gs_encode.sh

gsea_m="mean"
comb_0="c2.cp"
comb_1="c4_cancer_c2.cp"
script="gs_encode_single.py"
script="gs_encode.py"
training=$1 #T

module load anaconda

for target in docetaxel cisplatin bortezomib_a bortezomib_b #cisplatin #bortezomib_a bortezomib_b #docetaxel 
do
	for gsea in $comb_1 c4_cancer c4 #$comb_0 $comb_1
	do
		for norm in none sample gene batch
		do
			for g_filter in {20..200..20} #20 #30 40 50 100 150 200 #10 20 30 40 50 100 150 200 #50 100 200 500 # 20 40 50
			do
				for error in nrmse
				do

					if [ "$script" == "gs_encode.py" ]
						then 
						script_name="/tigress/zamalloa/GIT/gs_encode.py"
						file_name="${target} ${gsea} ${norm} ${gsea_m} ${g_filter} ${error}"

						export script_name file_name

						sbatch /tigress/zamalloa/GIT/gs_encode.cmd

					elif [ "$script" == "gs_encode_single.py" ]
						then
						if [ "$training" == "T" ]
							then
							script_name="/tigress/zamalloa/GIT/gs_encode_single.py"
							python GIT/gs_single_prep.py $target $gsea $norm $gsea_m $g_filter $error

							gsea_call="GSEA_FILES/SINGLE_GS/${target}_gee_expFALSE_geeprocexpTRUE_geetargetTRUE_geetargetprocexpTRUE_gsea${gsea}_genenorm${norm}_featmethod${gsea_m}_gfilter_${g_filter}_call"
							
							for g in `cat $gsea_call`
							do
								script_name="/tigress/zamalloa/GIT/gs_encode_single.py"
								file_name="${target} ${gsea} ${norm} ${gsea_m} ${g_filter} ${error} ${g}"
								
								export script_name file_name
								sbatch /tigress/zamalloa/GIT/gs_encode.cmd
							done

						elif [ "$training" == "F" ]
							then
							python GIT/gs_single_post.py $target $gsea $norm $gsea_m $g_filter $error
						fi
					fi

					echo "Done sending Autoencoder job"
				done
			done
		done
	done
done