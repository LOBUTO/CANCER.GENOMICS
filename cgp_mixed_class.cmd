#!/bin/bash
# parallel job using 1 processor and runs for 4:00 hours:
#SBATCH -t 65:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:1
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-user=zamalloa@princeton.edu

module load cudatoolkit/8.0
module load cudann/cuda-8.0
module load python

THEANO_FLAGS='mode=FAST_RUN,allow_gc=True,linker=c,device=gpu,lib.cnmem=1,floatX=float32,nvcc.fastmath=True' python ${script_name} ${file_name}
