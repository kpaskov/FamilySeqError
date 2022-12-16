# FamilySeqError

## Purpose
This project contains code for estimating error rates in sequencing data using family structure, as published in

Paskov K, Jung JY, Chrisman B, Stockham NT, Washington P, Varma M, Sun MW, Wall DP. Estimating sequencing error rates using families. BioData Mining. 2021 Dec;14(1):1-9.

## Input and output
This code starts with VCF files (with accompanying .tbi files), split by chromosome. It produces a JSON file containing highly granular sequencing error rates for each individual along with several other error metrics. Error rates are given for the probability of observing a variant call $v \in \\{0/0, 0/1, 1/1, ./.\\}$ given that the true genotype is $g \in \\{0/0, 0/1, 1/1\\}$. The JSON file contains

- $-log10(P[\text{obs}=v \mid \text{true\\_gen}=g])$
- $\text{lower\\_bound}[-log10(P[\text{obs}=v \mid \text{true\\_gen}=g])]$ (based on the number of observations in the data)
- $E[\text{obs}=v, \text{true\\_gen}=g]$
- $\text{precision}_g$
- $\text{recall}_g$
- $\text{FPR}_g$ (false positive rate)
- $\text{FNR}_g$ (false negative rate)
- $\text{F1}_g$ (F1-score)
- $\text{upper\\_bound}[\text{precision}_g]$ (based on the number of observations in the data)
- $\text{upper\\_bound}[\text{recall}_g]$ (based on the number of observations in the data)
- $\text{upper\\_bound}[\text{F1}_g]$ (based on the number of observations in the data)

Not all of these error metrics are available for parents. In particular, $P[\text{obs}=0/1 \mid \text{true\\_gen}=0/0]$ and $P[\text{obs}=0/1 \mid \text{true\\_gen}=1/1]$ can only be estimated for children (not parents) using this method.

## Instructions for running code
1. Start by getting your genomic data into numpy format using https://github.com/kpaskov/VCFtoNPZ. 

2. Pull family genotype counts.
A family genotype is a tuple of genotypes, representing the genotypes of a mother, father, and their child(ren), respectively, at a given site. The following code counts the number of times each family genotype occurs for each family on each chromosome.

`python parameter_estimation/pull_famgen_counts.py [data_dir] [ped_file] [chrom] [output_dir]`

3. Estimate sequencing error rates.
Now we can estimate error rates for each individual. Error rates are written to output_file in .json format.

`python parameter_estimation/estimate_parameters_per_individual.py [data_dir] [output_file]`
