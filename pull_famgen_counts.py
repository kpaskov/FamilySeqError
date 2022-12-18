import sys
from itertools import product
import os
import numpy as np
import scipy.sparse as sparse
import argparse
import json

parser = argparse.ArgumentParser(description='Pull family genotype counts.')
parser.add_argument('data_dir', type=str, help='Data directory of genotype data in .npy format produced using the VCFtoNPZ project.')
parser.add_argument('ped_file', type=str, help='Pedigree file (.ped).')
parser.add_argument('chrom', type=str, help='Chromosome.')
parser.add_argument('--count_type', type=str, default=None, help='Name of count type. Used to differentiate between counts in high-complexity vs low-complexity regions, for example.')
parser.add_argument('--include', type=str, default=None, nargs='+', help='Regions to include (.bed).')
parser.add_argument('--exclude', type=str, default=None, nargs='+', help='Regions to exclude (.bed).')
parser.add_argument('--use_pass', action='store_true', default=False, help='If flag is present, use apply pass filter to filter out variants that do not PASS. If flag is absent, ignore PASS filter.')
#parser.add_argument('--use_bases', action='store_true', default=False, help='Pull counts per base (ex. AA, AT) rather than per genotype (ex. 0/0, 0/1).')
args = parser.parse_args()

if not os.path.exists('%s/family_genotype_counts' % args.data_dir):
    os.makedirs('%s/family_genotype_counts' % args.data_dir)

if args.count_type is None:
    output_dir = '%s/family_genotype_counts'
else:
    output_dir = '%s/family_genotype_counts/%s' % (args.data_dir, args.count_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

if args.chrom == '23':
    args.chrom = 'X'
if args.chrom == '24':
    args.chrom = 'Y'
if args.chrom == '25':
    args.chrom = 'MT'

sample_file = '%s/genotypes/samples.json' % args.data_dir

obss = ['0/0', '0/1', '1/1', './.']

with open('%s/genotypes/info.json' % args.data_dir, 'r') as f:
    info = json.load(f)

info['ped_file'] = args.ped_file
info['include'] = args.include
info['exclude'] = args.exclude
info['use_pass'] = args.use_pass

with open('%s/info.json' % output_dir, 'w+') as f:
    json.dump(info, f)

def process_bedfile(bed_file):
    regions = []
    coverage = 0
    with open(bed_file, 'r') as f:
        for line in f:
            if '\t' in line.strip():
                pieces = line.strip().split('\t')
            else:
                pieces = line.strip().split(':')
                pieces = [pieces[0]] + pieces[1].strip().split('-')

            if pieces[0] == args.chrom or pieces[0] == 'chr%s' % args.chrom:
                regions.append(int(pieces[1]))
                regions.append(int(pieces[2])+1)
                coverage += (int(pieces[2])+1 - int(pieces[1]))
    return np.array(regions), coverage

if args.include is not None:
    include_regions = []
    for include_file in args.include:
        regions, coverage = process_bedfile(include_file)
        include_regions.append(regions)
        print('including %d bp (%s)' % (coverage, include_file))
else:
    include_regions = None

if args.exclude is not None:
    exclude_regions = []
    for exclude_file in args.exclude:
        regions, coverage = process_bedfile(exclude_file)
        exclude_regions.append(regions)
        print('excluding %d bp (%s)' % (coverage, exclude_file))
else:
    exclude_regions = None

# pull families with sequence data
with open(sample_file, 'r') as f:
    sample_id_to_index = dict([(sample_id, i) for i, sample_id in enumerate(json.load(f))])
# pull families from ped file
families = dict()
with open(args.ped_file, 'r') as f:	
    for line in f:
        pieces = line.strip().split('\t')
        if len(pieces) < 4:
            print('ped parsing error', line)
        else:
            fam_id = pieces[0]
            child_id, f_id, m_id = [x.replace('.', '_') for x in pieces[1:4]]

            if child_id in sample_id_to_index and f_id in sample_id_to_index and m_id in sample_id_to_index:
                if (fam_id, m_id, f_id) not in families:
                    families[(fam_id, m_id, f_id)] = [m_id, f_id]
                families[(fam_id, m_id, f_id)].append(child_id)
print('families %d' % len(families))


gen_files = sorted([f for f in os.listdir('%s/genotypes' % args.data_dir) if ('chr.%s.' % args.chrom) in f and 'gen.npz' in f], key=lambda x: int(x.split('.')[2]))
for gen_file in gen_files:
    batch_num = int(gen_file.split('.')[2])

    out_file = '%s/chr.%s.%d.famgen.counts.txt' % (output_dir, args.chrom, batch_num)
    print('saving to %s' % out_file)

    pos_data = np.load('%s/genotypes/chr.%s.%d.gen.coordinates.npy' % (args.data_dir, args.chrom, batch_num))
    if pos_data.shape[0]>0:
        is_snp = pos_data[:, 2].astype(bool)
        is_pass = pos_data[:, 3].astype(bool)

        is_ok = is_snp
        print('starting SNPs: %d' % np.sum(is_ok))

        if args.use_pass:
            is_ok = is_ok & is_pass
        print('surviving pass filter: %d' % np.sum(is_ok))
        
        if include_regions is not None:
            for regions, filename in zip(include_regions, args.include):
                insert_loc = np.searchsorted(regions, pos_data[:, 1])
                is_ok_include = np.remainder(insert_loc, 2)==1
                is_ok = is_ok & is_ok_include
                print('surviving include filter: %d (%s)' % (np.sum(is_ok), filename))

        if exclude_regions is not None:
            for regions, filename in zip(exclude_regions, args.exclude):
                insert_loc = np.searchsorted(regions, pos_data[:, 1])
                is_ok_exclude = np.remainder(insert_loc, 2)==0
                is_ok = is_ok & is_ok_exclude
                print('surviving exclude filter: %d (%s)' % (np.sum(is_ok), filename))

        # Pull data together
        A = sparse.load_npz('%s/genotypes/%s' % (args.data_dir, gen_file))

        # filter out snps
        A = A[:, is_ok]
        print('genotype matrix prepared', A.shape)
    else:
        A = np.zeros((0, 0))

    with open(out_file, 'w+') as f: 
        for famkey, inds in families.items():
            m = len(inds)
            genotype_to_counts = np.zeros((4,)*m, dtype=int)

            if A.shape[1]>0:    
                indices = [sample_id_to_index[ind] for ind in inds]
                family_genotypes = A[indices, :]
                    
                # remove positions where whole family is homref
                has_data = sorted(set(family_genotypes.nonzero()[1]))
                num_hom_ref = family_genotypes.shape[1] - len(has_data)

                family_genotypes = family_genotypes[:, has_data].A
                #print(famkey, family_genotypes.shape)

                # recode missing values
                family_genotypes[family_genotypes<0] = 3

                #print(A.shape, family_genotypes.shape)
                    
                # fill in genotype_to_counts
                unique_gens, counts = np.unique(family_genotypes, return_counts=True, axis=1)
                for g, c in zip(unique_gens.T, counts):
                    genotype_to_counts[tuple(g)] += c

                # add hom ref sites (that were previously removed)
                genotype_to_counts[(0,)*m] = num_hom_ref

            # write to file
            f.write('\t'.join([famkey[0], '.'.join(inds)] + \
                [str(genotype_to_counts[g]) for g in product([0, 1, 2, 3], repeat=m)]) + '\n')
            
