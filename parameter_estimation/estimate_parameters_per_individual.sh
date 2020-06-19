#!/bin/bash
#
#
#SBATCH --job-name=params
#SBATCH --output=logs/params.out
#SBATCH --error=logs/params.err
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 30:00:00
#SBATCH --mem=8G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

srun python3 parameter_estimation/estimate_parameters_per_individual.py ../DATA/spark/family_genotype_counts/quads params/spark_quads_params.json
