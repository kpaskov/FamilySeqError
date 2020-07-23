import numpy as np
from scipy.sparse import csc_matrix, save_npz, hstack
import time
import argparse
import gzip
from pysam import VariantFile, TabixFile
import json
import os
import itertools

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('vcf_file', type=str, help='VCF file to pull from.')
parser.add_argument('assembly', type=str, help='Human genome reference used.')
parser.add_argument('out_directory', type=str, help='Output directory.')
parser.add_argument('chrom', type=str, help='Chromosome of interest.')
parser.add_argument('--batch_size', type=int, default=-1, help='Restrict number of positions per file to batch_size.')
parser.add_argument('--batch_num', type=int, default=0, help='To be used along with batch_size to restrict positions per file. Will include positions >= batch_num*batch_size and <= (batch_num+1)*batch_size')
parser.add_argument('--maxsize', type=int, default=500000000, help='Amount of memory per block.')
parser.add_argument('--additional_vcf_files', type=str, nargs='+', help='Additional VCF files to pull data from.')
parser.add_argument('--id_mapper_file', type=str, default=None, help='File that maps old ids to new ones.')
parser.add_argument('--id_mapper_sep', type=str, default='\t', help='Separater to parse id_mapper_file.')
parser.add_argument('--old_id_index', type=int, default=0, help='Index of old_id in id_mapper_file.')
parser.add_argument('--new_id_index', type=int, default=1, help='Index of new_id in id_mapper_file.')
args = parser.parse_args()

t0 = time.time()

chrom_int = 23 if args.chrom == 'X' else 24 if args.chrom == 'Y' else 25 if args.chrom == 'MT' else int(args.chrom)

gen_mapping = {'./.': -1, '0/0': 0, '0|0': 0, '0/1': 1, '0|1': 1, '1/0': 1, '1|0': 1, '1/1': 2, '1|1': 2}

def process_header(vcf):
    sample_ids = [x.replace('.', '_') for x in vcf.header.samples]

    if args.id_mapper_file is not None:
        old_id_to_new_id = dict()
        with open(args.id_mapper_file, 'r') as f:
            for line in f:
                pieces = line.strip().split(args.id_mapper_sep)
                if len(pieces)>args.old_id_index and len(pieces)>args.new_id_index:
                    old_id_to_new_id[pieces[args.old_id_index]] = pieces[args.new_id_index]
        sample_ids = [old_id_to_new_id[x] for x in sample_ids]

    sample_file = '%s/samples.json' % args.out_directory
    if os.path.isfile(sample_file):
        with open(sample_file, 'r') as f:
            stored_sample_ids = json.load(f)
            assert sample_ids == stored_sample_ids
    else:
        with open(sample_file, 'w+') as f:
            json.dump(sample_ids, f)

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
    np.save('%s/chr.%s.%d.gen.coordinates' % (args.out_directory, args.chrom, args.batch_num), np.asarray(np.asarray(chrom_coord, dtype=int), dtype=int))
    print('Completed in ', time.time()-t0, 'sec')

with open('%s/info.json' % args.out_directory, 'w+') as f:
    json.dump({'assembly': args.assembly, 'batch_size': args.batch_size, 'vcf_directory': '/'.join(args.vcf_file.split('/')[:-1])}, f)

vcf = VariantFile(args.vcf_file)
sample_ids, contigs = process_header(vcf)

if len(args.additional_vcf_files)>0:
    for vcf_file in args.additional_vcf_files:
        if os.path.isfile(vcf_file):
            new_vcf = VariantFile(vcf_file)
            new_sample_ids, _ = process_header(new_vcf)
            assert sample_ids == new_sample_ids
        else:
            print(vcf_file, 'does not exist')

contig = None
if args.chrom in contigs:
    contig = contigs[args.chrom]
elif 'chr%s' % args.chrom in contigs:
    contig = contigs['chr%s' % args.chrom]
else:
    raise Exception('Trouble finding contig', args.chrom, 'in', contig_names)
print('Chrom length', contig.length)

vcfs = [TabixFile(args.vcf_file, parser=None)]
for vcf_file in args.additional_vcf_files:
    if os.path.isfile(vcf_file):
        vcfs.append(TabixFile(vcf_file, parser=None))
    else:
        print(vcf_file, 'does not exist')

if args.batch_size != -1:
    start_pos, end_pos = args.batch_num*args.batch_size, (args.batch_num+1)*args.batch_size
    print('Interval', start_pos, end_pos)
    if start_pos < contig.length:
        process_body(itertools.chain(*[vcf.fetch(reference=contig.name, start=start_pos, end=end_pos) for vcf in vcfs]), sample_ids)
    else:
        print('Interval (%d-%d) is longer than chromosome (length=%d).' % (start_pos, end_pos, contig.length))
else:
    process_body(itertools.chain(*[vcf.fetch(reference=contig.name) for vcf in vcf]), sample_ids)




