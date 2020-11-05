#!/bin/bash
#
#
#SBATCH --job-name=pull
#SBATCH --output=logs/pull.out
#SBATCH --error=logs/pull.err
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 10:00:00
#SBATCH --mem=64G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
module load biology
module load py-pysam/0.15.3_py36

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID is " $SLURM_ARRAY_TASK_ID

srun python3 preprocessing/pull_gen_data.py /scratch/PI/dpwall/DATA/iHART/vcf/mssng/mssng_db6.chr$1.recal.vcf.gz ../DATA/mssng/genotypes $1 --batch_size 10000000 --batch_num $2 --maxsize 5000000000

#srun python3 preprocessing/pull_gen_data.py /oak/stanford/groups/dpwall/simons_spark/snp/SPARK.30K.array_genotype.20190818.phaseable.passing.vcf ../DATA/spark/genotypes X

#srun python3 preprocessing/pull_gen_data.py /oak/stanford/groups/dpwall/simons_spark/snp/SPARK.30K.array_genotype.20190818.phaseable.passing.vcf ../DATA/spark/genotypes Y
