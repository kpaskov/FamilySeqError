#!/bin/bash
#
#
#SBATCH --job-name=pass
#SBATCH --output=logs/pass.out
#SBATCH --error=logs/pass.err
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 1:00:00
#SBATCH --mem=8G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36


srun python3 preprocessing/pull_pass.py ../DATA/spark/genotypes --pass_from_gen ../DATA/spark/spark.ped
