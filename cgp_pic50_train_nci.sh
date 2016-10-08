#cgp_pic50_train_nci.sh
#Master script to build pre training files to target specific NSC compound

Rscript /tigress/zamalloa/GIT/cgp_pic50_train_nci.R $1

module load cudatoolkit
module load python
python /tigress/zamalloa/GIT/cgp_pic50_train.py $1 #Same function as before sufficient?

echo "Done"
