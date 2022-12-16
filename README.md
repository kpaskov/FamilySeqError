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

The code adds to a directory structure created by the https://github.com/kpaskov/VCFtoNPZ project.

```
data_dir
- genotypes
- family_genotype_counts
- - info.json
- - chr.1.0.famgen.counts.txt
- - chr.1.1.famgen.counts.txt
...
- sequencing_error_rates
- - errors.json
```

The `info.json` file contains metadata including the reference assembly (GRch37 or GRch38) used to produce the variant calls, the vcf_directory used to generate the data, and if relevant the batch_size used to generate the data. It also includes the path to the ped_file used to generate the counts, as well as the variant filters applied.

The `chr.[chrom].[batch_num].famgen.counts.txt` files contain counts of the number of times each family genotype combination occurs within the corresponding segment of chromosome `[chrom]`, with variant filters applied. The `[data_dir]/family_genotype_counts` directory may contain subdirectories corresponding to different variant filters (see `--count_type` option below).

The `errors.json` file contains the sequencing error rates estimated by our method. There may be multiple files in the `[data_dir]/sequencing_error_rates` directory corresponding to different variant filters. For example, our method can be used to estimate sequencing error rates in low-complexity vs high-complexity regions (again, see `--count_type` option below).

## Instructions for running code (autosomal error rates)
1. Start by getting your genomic data into numpy format using https://github.com/kpaskov/VCFtoNPZ. 

2. Pull family genotype counts.
A family genotype is a tuple of genotypes, representing the genotypes of a mother, father, and their child(ren), respectively, at a given site. The following code counts the number of times each family genotype occurs for each family on each segment of each chromosome. It produces a series of files `[output_dir]/chr.[chrom].[batch_num].famgen.counts.txt` which contain a line for each family in the dataset, representing the number of times each family genotype occurs for that family within the corresponding segment of chromosome `[chrom]`.

`python pull_famgen_counts.py [data_dir] [ped_file] [chrom]`

The script has options 
- `--use_pass` which uses the PASS flag (from the VCF file) to filter variants. Only variants that PASS are counted.
- `--include [bed_file]` which uses a BED file to filter variants. Only variants within the intervals listed in the BED file are counted.
- `--exlude [bed_file]` which uses a BED file to filter variants. Only variants outside the intervals listed in the BED file are counted.
- `--exlude [bed_file]` which uses a BED file to filter variants. Only variants outside the intervals listed in the BED file are counted.
- `--count_type [output_dir]` which creates a subdirectory `[data_dir]/family_genotype_counts/[output_dir]` and stores the `*.famgen.counts.txt` files there. This option is useful when creating different types of family genotype counts from the same dataset, for example counts in low-complexity and high-complexity regions.

3. Estimate sequencing error rates.
Now we can estimate error rates for each individual. Error rates are written to `[output_file]` in .json format.

`python estimate_sequencing_error_rates.py [data_dir]`

The script has options 
- `--is_ngs` to be used when working with next-gen sequencing data.
- `--count_type [input_dir]` which uses the family genotype counts stored in `[data_dir]/family_genotype_counts/[input_dir]` to estimate sequencing error rates. This option is useful when creating different types of family genotype counts from the same dataset, for example counts in low-complexity and high-complexity regions.

## Instructions for estimating error rates on the X-chromosome

