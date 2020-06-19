#!/bin/bash
#
#
#SBATCH --job-name=pull
#SBATCH --output=logs/pull.out
#SBATCH --error=logs/pull.err
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 10:00:00
#SBATCH --mem=8G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID is " $SLURM_ARRAY_TASK_ID

#srun python3 preprocessing/pull_gen_data.py /scratch/PI/dpwall/DATA/iHART/vcf/mssng/mssng_db6.chr$SLURM_ARRAY_TASK_ID.recal.vcf.gz ../DATA/mssng/genotypes $SLURM_ARRAY_TASK_ID --batch_size 1000000 --batch_num $1

srun python3 preprocessing/pull_gen_data.py /oak/stanford/groups/dpwall/simons_spark/snp/SPARK.30K.array_genotype.20190818.phaseable.passing.vcf ../DATA/spark/genotypes X

srun python3 preprocessing/pull_gen_data.py /oak/stanford/groups/dpwall/simons_spark/snp/SPARK.30K.array_genotype.20190818.phaseable.passing.vcf ../DATA/spark/genotypes Y
