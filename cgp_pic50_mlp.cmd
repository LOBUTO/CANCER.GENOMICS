#!/bin/bash
# parallel job using 1 processor and runs for 2:00 hours:
#SBATCH -t 5:00:00
#SBATCH --mem-per-cpu=15000
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:1
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-user=zamalloa@princeton.edu

module load cudatoolkit
module load python

THEANO_FLAGS='device=gpu' python GIT/cgp_pic50_mlp.py 50_50_50 Erlotinib
#THEANO_FLAGS='device=gpu' python GIT/cgp_test.py
#THEANO_FLAGS='device=gpu' python GIT/cgp_pic50_testing.py CGP_FILES/nci60_cgp_cor.csv CGP_FILES/CGP_RESULTS/Erlotinib_cgp_pIC50.200.200.pkl CGP_FILES/CGP_TRAIN_TABLES/STD_SCALER.Erlotinib.pIC50.pkl CGP_FILES/nci60_cgpfeat_cancertable.csv Erlotinib 0.7
