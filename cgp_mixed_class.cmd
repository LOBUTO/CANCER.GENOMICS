#!/bin/bash
# parallel job using 1 processor and runs for 4:00 hours:
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu=6000
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

THEANO_FLAGS='device=gpu' python GIT/cgp_mixed_class.py {file_name}
