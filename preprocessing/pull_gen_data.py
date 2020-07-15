import numpy as np
from scipy.sparse import csc_matrix, save_npz
import time
import argparse
import gzip
from pysam import VariantFile, TabixFile
import json

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('vcf_file', type=str, help='VCF file to pull from.')
parser.add_argument('assembly', type=str, help='Human genome reference used.')
parser.add_argument('out_directory', type=str, help='Output directory.')
parser.add_argument('chrom', type=str, help='Chromosome of interest.')
parser.add_argument('--batch_size', type=int, default=None, help='Restrict number of positions per file to batch_size.')
parser.add_argument('--batch_num', type=int, default=0, help='To be used along with batch_size to restrict positions per file. Will include positions >= batch_num*batch_size and <= (batch_num+1)*batch_size')
parser.add_argument('--maxsize', type=int, default=500000000, help='Amount of memory per block.')
args = parser.parse_args()

t0 = time.time()

chrom_int = 23 if args.chrom == 'X' else 24 if args.chrom == 'Y' else 25 if args.chrom == 'MT' else int(args.chrom)

gen_mapping = {'./.': -1, '0/0': 0, '0|0': 0, '0/1': 1, '0|1': 1, '1/0': 1, '1|0': 1, '1/1': 2, '1|1': 2}

def process_header(vcf):
    sample_ids = [x.replace('.', '_') for x in vcf.header.samples]

    if args.batch_num == 0:
        with open('%s/chr.%s.gen.samples.txt' % (args.out_directory, args.chrom), 'w+') as sample_f:
            # Pull sample_ids and write to file
            sample_f.write('\n'.join(sample_ids))
            print('Num individuals with genomic data', len(sample_ids))
    return sample_ids, vcf.header.contigs

def process_body(records, sample_ids):

    data, indices, indptr, index = np.zeros((args.maxsize,), dtype=np.int8), np.zeros((args.maxsize,), dtype=int), [0], 0
    chrom_coord = []

    with gzip.open('%s/chr.%s.%d.gen.variants.txt.gz' % (args.out_directory, args.chrom, args.batch_num), 'wt') as variant_f:
        for line in records:
            pieces = line.strip().split('\t')
            fmt = pieces[8].strip().split(':')

            # Write variant to file
            variant_f.write('\t'.join(pieces[:9]) + '\n')

            # pull chrom_coord information
            pos, _, ref, alt = pieces[1:5]
            is_biallelic_snp = 1 if len(ref) == 1 and len(alt) == 1 and ref != '.' and alt != '.' else 0
            is_pass = pieces[6] == 'PASS'
            chrom_coord.append((chrom_int, int(pos), is_biallelic_snp, is_pass))

            # pull genotypes
            gen_index = fmt.index('GT')
            for i, piece in enumerate(pieces[9:]):
                segment = piece.split(':', maxsplit=gen_index+1)
                gt = gen_mapping.get(segment[gen_index], -1) # For now we mark multi-base loci as unknown

                if gt != 0:
                    indices[index] = i
                    data[index] = gt
                    index += 1
            indptr.append(index)

    gen = csc_matrix((data[:index], indices[:index], indptr), shape=(len(sample_ids), len(indptr)-1), dtype=np.int8)
        
    # Save to file
    save_npz('%s/chr.%s.%d.gen' % (args.out_directory, args.chrom, args.batch_num), gen)
    np.save('%s/chr.%s.%d.gen.coordinates' % (args.out_directory, args.chrom, args.batch_num), np.asarray(chrom_coord, dtype=int))
    print('Completed in ', time.time()-t0, 'sec')

with open('%s/info.json' % args.out_directory, 'w+') as f:
    json.dump({'assembly': args.assembly, 'batch_size': args.batch_size, 'vcf_directory': '/'.join(args.vcf_file.split('/')[:-1])}, f)

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

vcf = TabixFile(args.vcf_file, parser=None)
if args.batch_size is not None:
    start_pos, end_pos = args.batch_num*args.batch_size, (args.batch_num+1)*args.batch_size
    print('Interval', start_pos, end_pos)
    if start_pos < contig.length:
        process_body(vcf.fetch(reference=contig.name, start=start_pos, end=end_pos), sample_ids)
    else:
        print('Interval (%d-%d) is longer than chromosome (length=%d).' % (start_pos, end_pos, contig.length))
else:
    process_body(vcf.fetch(reference=contig.name), sample_ids)

