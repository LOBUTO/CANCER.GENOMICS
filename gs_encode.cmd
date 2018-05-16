#!/bin/bash
# parallel job using 1 processor and runs for 4:00 hours:
#SBATCH -t 3:00:00
#SBATCH --mem-per-cpu=32000
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-socket=1
#SBATCH --gres=gpu:2
#SBATCH --exclude=tiger-i20g7,tiger-i21g3,tiger-i22g3,tiger-i21g1
# sends mail when process begins, and
# when it ends. Make sure you define your email
# address.
#SBATCH --mail-user=zamalloa@princeton.edu

module load anaconda
module load cudatoolkit/8.0
module load cudann/cuda-8.0

python ${script_name} ${file_name}