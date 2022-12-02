# FamilySeqError

## Purpose
This project contains code for estimating error rates in sequencing data using family structure.

## Input and output

## Instructions for running code
1. Start by getting your genomic data into numpy format using [(https://github.com/kpaskov/VCFtoNPZ)](https://github.com/kpaskov/VCFtoNPZ). 

2. Pull family genotype counts.
A family genotype is a tuple of genotypes, representing the genotypes of a mother, father, and their child(ren), respectively, at a given site. The following code counts the number of times each family genotype occurs for each family on each chromosome.

`python parameter_estimation/pull_famgen_counts.py [data_dir] [ped_file] [chrom] [output_dir]`

3. Estimate sequencing error rates.
Now we can estimate error rates for each individual. Error rates are written to output_file in .json format.

`python parameter_estimation/estimate_parameters_per_individual.py [data_dir] [output_file]`
