#!/bin/bash
#cgp_pic50_train.sh


Rscript /tigress/zamalloa/GIT/cgp_pic50_train.R $1

module load cudatoolkit
module load python
python /tigress/zamalloa/GIT/cgp_pic50_train.py $1

echo "Done"
