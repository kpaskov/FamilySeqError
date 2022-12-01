This project contains code for estimating error rates in sequencing data using family structure.

1. Preprocessing.
Start by getting your genomic data into numpy format. If your data is currently in VCF format, split by chromosome, this can be done by running

python preprocessing/pull_gen_data.py [vcf_file] [data_dir] [chrom]

If your vcf files don't have filters applied (for example no variant is PASS) or you'd like to apply a different type of filter, use preprocessing/pull_pass.py

2. Pull family genotype counts.
A family genotype is a tuple of genotypes, representing the genotypes of a mother, father, and their child(ren), respectively, at a given site. The following code counts the number of times each family genotype occurs for each family on each chromosome.

python parameter_estimation/pull_famgen_counts.py [data_dir] [ped_file] [chrom] [output_dir]

3. Estimate sequencing error rates.
Now we can estimate error rates for each individual. Error rates are written to output_file in .json format.

python parameter_estimation/estimate_parameters_per_individual.py [data_dir] [output_file]

fixed [('AU035806', 'mat'), ('AU035806', 'pat'), ('AU035810', 'mat'), ('AU035810', 'pat'), ('AU035812', 'mat'), ('AU035812', 'pat')]
