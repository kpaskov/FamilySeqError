#!/bin/bash
#
#
#SBATCH --job-name=pull
#SBATCH --output=logs/pull.out
#SBATCH --error=logs/pull.err
#SBATCH --array=0-0
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 3:00:00
#SBATCH --mem=16G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
module load biology
module load py-pysam/0.15.3_py36

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID is " $SLURM_ARRAY_TASK_ID

srun python3 preprocessing/pull_gen_data.py /scratch/groups/dpwall/DATA/SSC/SPARK/joint.vcf/merged_WES_v0.vcf.gz 

#srun python3 preprocessing/pull_gen_data.py /scratch/PI/dpwall/DATA/iHART/vcf/mssng/mssng_db6.chr$SLURM_ARRAY_TASK_ID.recal.vcf.gz 38 ../DATA/mssng/genotypes $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 1300000000

#srun python3 preprocessing/pull_gen_data.py /oak/stanford/groups/dpwall/simons_spark/snp/SPARK.30K.array_genotype.20190818.phaseable.passing.vcf 38 ../DATA/spark/genotypes $SLURM_ARRAY_TASK_ID

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/v34/$SLURM_ARRAY_TASK_ID.reheader.vcf.gz 37 ../DATA/ihart/genotypes $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc/ssc.$SLURM_ARRAY_TASK_ID.vcf.gz 37 ../DATA/ssc/genotypes $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-1.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-1 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-2.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-2 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-3.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-3 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-4.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-4 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-5.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-5 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase1-7.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase1-7 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase2.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase2 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase2_B01.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase2_B01 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase2_Replacements.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase2_Replacements $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase3_1.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase3_1 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase3_1_B02.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase3_1_B02 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase3_2.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase3_2 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/phase4.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/phase4 $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py $OAK/ihart/vcf/ssc.hg38/pilot.all.vcf.gz 38 ../DATA/ssc.hg38/genotypes/pilot $SLURM_ARRAY_TASK_ID --batch_size 10000000 --batch_num $1 --maxsize 500000000

#srun python3 preprocessing/pull_gen_data.py /scratch/PI/dpwall/DATA/iHART/vcf/ms2.hg38/chr$1.01.vcf.gz 38 ../DATA/ihart.ms2/genotypes $1 --batch_size 10000000 --batch_num $SLURM_ARRAY_TASK_ID --maxsize 1300000000 --additional_vcf_files /scratch/PI/dpwall/DATA/iHART/vcf/ms2.hg38/chr$1.02.vcf.gz /scratch/PI/dpwall/DATA/iHART/vcf/ms2.hg38/chr$1.03.vcf.gz /scratch/PI/dpwall/DATA/iHART/vcf/ms2.hg38/chr$1.04.vcf.gz /scratch/PI/dpwall/DATA/iHART/vcf/ms2.hg38/chr$1.05.vcf.gz --id_mapper_file data/ihart_id_map.txt
