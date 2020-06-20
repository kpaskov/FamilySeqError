from collections import defaultdict
import numpy as np
from scipy.sparse import csc_matrix, save_npz
import time
from itertools import product, islice
import sys
import argparse
import gzip
from pysam import VariantFile
from collections.abc import Iterable

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('vcf_file', type=str, help='VCF file to pull from.')
parser.add_argument('out_directory', type=str, help='Output directory.')
parser.add_argument('chrom', type=str, help='Chromosome of interest.')
parser.add_argument('--batch_size', type=int, default=None, help='Restrict number of positions per file to batch_size.')
parser.add_argument('--batch_num', type=int, default=0, help='To be used along with batch_size to restrict positions per file. Will include positions >= batch_num*batch_size and <= (batch_num+1)*batch_size')
parser.add_argument('--maxsize', type=int, default=500000000, help='Amount of memory per block.')
args = parser.parse_args()

t0 = time.time()

chrom_int = 23 if args.chrom == 'X' else 24 if args.chrom == 'Y' else 25 if args.chrom == 'MT' else int(args.chrom)

gen_mapping = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 2}

def process_header(vcf):
    sample_ids = vcf.header.samples

    if args.batch_num == 0:
        with open('%s/chr.%s.gen.samples.txt' % (args.out_directory, args.chrom), 'w+') as sample_f:
            # Pull sample_ids and write to file
            sample_f.write('\n'.join(sample_ids))
            print('Num individuals with genomic data', len(sample_ids))
    return sample_ids, vcf.header.contigs

def process_body(records, sample_ids):

    data, indices, indptr, index = np.zeros((args.maxsize,), dtype=np.int8), np.zeros((args.maxsize,), dtype=int), [0], 0

    with gzip.open('%s/chr.%s.%d.gen.variants.txt.gz' % (args.out_directory, args.chrom, args.batch_num), 'wt') as variant_f:
        num_lines_in_file = 0
        chrom_coord = []
        for record in records:

            # Write variant to file
            variant_f.write('\t'.join(str(record).split(maxsplit=9)[:9]) + '\n')

            is_biallelic_snp = 1 if len(record.ref) == 1 and len(record.alts) == 1 and len(record.alts[0]) == 1 and record.ref != '.' and record.alts[0] != '.' else 0
            is_pass = 'PASS' in record.filter.keys()
            chrom_coord.append((chrom_int, int(record.pos), is_biallelic_snp, is_pass))

            # Pull out genotypes
            for i, sample_id in enumerate(sample_ids):
                gt = gen_mapping.get(record.samples[sample_id]['GT'], -1)
                if gt != 0:
                    indices[index] = i
                    data[index] = gt
                    index += 1
            indptr.append(index)
            num_lines_in_file += 1


    gen = csc_matrix((data[:index], indices[:index], indptr), shape=(len(sample_ids), num_lines_in_file), dtype=np.int8)
        
    # Save to file
    save_npz('%s/chr.%s.%d.gen' % (args.out_directory, args.chrom, args.batch_num), gen)
    np.save('%s/chr.%s.%d.gen.coordinates' % (args.out_directory, args.chrom, args.batch_num), np.asarray(chrom_coord, dtype=int))
    print('Completed in ', time.time()-t0, 'sec')

vcf = VariantFile(args.vcf_file)
sample_ids, contigs = process_header(vcf)

contig = None
if args.chrom in contigs:
    contig = contigs[args.chrom]
elif 'chr%s' % args.chrom in contigs:
    contig = contigs['chr%s' % args.chrom]
else:
    raise Exception('Trouble finding contig', args.chrom, 'in', contig_names)

print('Chrom length', contig.length)
if args.batch_size is not None:
    start_pos, end_pos = args.batch_num*args.batch_size, (args.batch_num+1)*args.batch_size
    print('Interval', start_pos, end_pos)
    if start_pos < contig.length:
        process_body(vcf.fetch(contig=contig.name, start=start_pos, stop=end_pos), sample_ids)
    else:
        print('Interval (%d-%d) is longer than chromosome (length=%d).' % (start_pos, end_pos, contig.length))
else:
    process_body(vcf.fetch(contig=contig.name), sample_ids)

