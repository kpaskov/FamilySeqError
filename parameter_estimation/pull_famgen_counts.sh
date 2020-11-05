#!/bin/bash
#
#
#SBATCH --job-name=famgen
#SBATCH --output=logs/famgen_%a.out
#SBATCH --error=logs/famgen_%a.err
#SBATCH --array=1-22
#SBATCH -p dpwall
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/FamilySeqError
#SBATCH -t 5:00:00
#SBATCH --mem=16G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID is " $SLURM_ARRAY_TASK_ID

# ------------------------ For parameter estimation paper -----------------------

# ihart
#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart.ms2/genotypes ../DATA/ihart.ms2/ihart.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart.ms2/family_genotype_counts

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart.ms2/genotypes ../DATA/ihart.ms2/ihart.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart.ms2/family_genotype_counts/HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart.ms2/genotypes ../DATA/ihart.ms2/ihart.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart.ms2/family_genotype_counts/LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

# ihart (GATK 3.2)
# srun python3 parameter_estimation/pull_famgen_counts.py split_gen_ihart3.2 data/v34.vcf.ped $SLURM_ARRAY_TASK_ID split_gen_ihart3.2

# ihart identicals
# srun python3 parameter_estimation/pull_famgen_counts.py split_gen_ihart data/ihart_identicals.ped $SLURM_ARRAY_TASK_ID split_gen_ihart_identicals

# spark exome EX
# srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX --include data/VCRome_2_1_hg38_primary_targets_liftover_ordered.bed 
#srun python3 parameter_estimation/pull_famgen_counts_vcf.py /oak/stanford/groups/dpwall/ihart/vcf/v34/$SLURM_ARRAY_TASK_ID.reheader.vcf.gz data/v34.vcf.ped $SLURM_ARRAY_TASK_ID counts_ihart --family_sizes 3 4 5 6 7 --depth_bins 20 30

# spark exome EX25
# srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX25 --include data/VCRome_2_1_hg38_primary_targets_liftover_0_25.bed

# spark exome EX50
# srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX50 --include data/VCRome_2_1_hg38_primary_targets_liftover_25_50.bed

# spark exome EX75
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX75 --include data/VCRome_2_1_hg38_primary_targets_liftover_50_75.bed

# spark exome EX1000
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX1000 --include data/VCRome_2_1_hg38_primary_targets_liftover_75_1000.bed

# spark exome identicals
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark_exome data/spark_exome_identicals.ped $SLURM_ARRAY_TASK_ID split_gen_spark_exome_EX_identicals --include data/VCRome_2_1_hg38_primary_targets_liftover_ordered.bed

# spark
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark data/spark.ped $SLURM_ARRAY_TASK_ID split_gen_spark

# spark identicals
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_spark data/spark_identicals.ped $SLURM_ARRAY_TASK_ID split_gen_spark_identicals

# ------------------------ For sibpair paper ------------------------
#srun python3 parameter_estimation/pull_famgen_counts.py ../PhasingFamilies/split_gen_ihart ../PhasingFamilies/data/v34.vcf.ped.quads.ped $SLURM_ARRAY_TASK_ID ../PhasingFamilies/split_gen_ihart_quads

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/spark/genotypes ../DATA/spark/spark.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/spark/family_genotype_counts/quads

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/spark/genotypes ../DATA/spark/spark.ped $SLURM_ARRAY_TASK_ID ../DATA/spark/family_genotype_counts

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/mssng/genotypes ../DATA/mssng/mssng_db6.vcf.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/mssng/family_genotype_counts/quads/hcr --exclude ../DATA/lcr_regions/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/mssng/genotypes ../DATA/mssng/mssng_db6.vcf.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/mssng/family_genotype_counts/quads/lcr --include ../DATA/lcr_regions/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ancestry/genotypes ../DATA/ancestry/ancestryDNA.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/ancestry/family_genotype_counts/quads

# ------------------------ For deletions paper -----------------------    

# platinum
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_platinum data/platinum.ped $SLURM_ARRAY_TASK_ID split_gen_platinum

# SSC
#srun python3 parameter_estimation/pull_famgen_counts.py split_gen_ssc data/ssc.ped $SLURM_ARRAY_TASK_ID split_gen_ssc

# ------------------------ For crossover paper ------------------------

# ihart
#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart/genotypes ../DATA/ihart/ihart.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart/family_genotype_counts/quads/LCR --include data/btu356-suppl_data/btu356_LCR-hs37d5.bed/btu356_LCR-hs37d5.bed 

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart/genotypes ../DATA/ihart/ihart.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart/family_genotype_counts/quads/HCR --exclude data/btu356-suppl_data/btu356_LCR-hs37d5.bed/btu356_LCR-hs37d5.bed    

# mssng
#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/mssng/genotypes ../DATA/mssng/mssng.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/mssng/family_genotype_counts/LCR_quads --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/mssng/genotypes ../DATA/mssng/mssng.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/mssng/family_genotype_counts/HCR_quads --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

# ssc
#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-1 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-1_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-1 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-1_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-2_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed                 

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-2_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-3 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-3_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-3 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-3_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-4 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-4_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-4 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-4_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-5 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-5_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-5 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-5_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-7 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-7_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase1-7 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase1-7_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2_B01 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_B01_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2_B01 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_B01_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2_Replacements ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_Replacements_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase2_Replacements ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase2_Replacements_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_1 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_1_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_1 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_1_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_1_B02 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_1_B02_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_1_B02 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_1_B02_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_2_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase3_2 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase3_2_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase4 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase4_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/phase4 ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/phase4_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/pilot ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/pilot_LCR --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ssc.hg38/genotypes/pilot ../DATA/ssc.hg38/ssc.ped $SLURM_ARRAY_TASK_ID ../DATA/ssc.hg38/family_genotype_counts/pilot_HCR --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart.ms2/genotypes ../DATA/ihart.ms2/ihart.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart.ms2/family_genotype_counts/HCR_quads --exclude data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed 

#srun python3 parameter_estimation/pull_famgen_counts.py ../DATA/ihart.ms2/genotypes ../DATA/ihart.ms2/ihart.ped.quads.ped $SLURM_ARRAY_TASK_ID ../DATA/ihart.ms2/family_genotype_counts/LCR_quads --include data/btu356-suppl_data/btu356_LCR-hs38.bed/btu356_LCR-hs38.bed
