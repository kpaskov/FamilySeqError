import sys
from itertools import product
from os import listdir
import numpy as np
import scipy.sparse as sparse
import argparse
import json

parser = argparse.ArgumentParser(description='Pull family genotype counts.')
parser.add_argument('data_dir', type=str, help='Directory of genotype data in .npy format.')
parser.add_argument('ped_file', type=str, help='Pedigree file (.ped).')
parser.add_argument('chrom', type=str, help='Chromosome.')
parser.add_argument('out_dir', type=str, help='Directory to write counts.')
parser.add_argument('--include', type=str, default=None, help='Regions to include (.bed).')
parser.add_argument('--exclude', type=str, default=None, help='Regions to exclude (.bed).')
#parser.add_argument('--use_bases', action='store_true', default=False, help='Pull counts per base (ex. AA, AT) rather than per genotype (ex. 0/0, 0/1).')
args = parser.parse_args()


if args.chrom == '23':
    args.chrom = 'X'
if args.chrom == '24':
    args.chrom = 'Y'
if args.chrom == '25':
    args.chrom = 'MT'

sample_file = '%s/samples.json' % args.data_dir

obss = ['0/0', '0/1', '1/1', './.']

with open('%s/info.json' % args.data_dir, 'r') as f:
    info = json.load(f)

with open('%s/info.json' % args.out_dir, 'w+') as f:
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

include_regions, include_coverage = (None, 0) if args.include is None else process_bedfile(args.include)
exclude_regions, exclude_coverage = (None, 0) if args.exclude is None else process_bedfile(args.exclude)

print('including %s bp' % ('all' if include_regions is None else str(include_coverage)))
print('excluding %s bp' % ('no' if exclude_regions is None else str(exclude_coverage)))

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

            if child_id in sample_id_to_index and (f_id in sample_id_to_index and m_id not in sample_id_to_index):
                if (fam_id, m_id, f_id) not in families:
                    families[(fam_id, m_id, f_id)] = [f_id]
                families[(fam_id, m_id, f_id)].append(child_id)
            elif child_id in sample_id_to_index and (m_id in sample_id_to_index and f_id not in sample_id_to_index):
                if (fam_id, m_id, f_id) not in families:
                    families[(fam_id, m_id, f_id)] = [m_id]
                families[(fam_id, m_id, f_id)].append(child_id)

families = dict([(k, v) for k, v in families.items() if len(v)==3])
print('families %d' % len(families))


gen_files = sorted([f for f in listdir(args.data_dir) if ('chr.%s.' % args.chrom) in f and 'gen.npz' in f], key=lambda x: int(x.split('.')[2]))
for gen_file in gen_files:
    batch_num = int(gen_file.split('.')[2])

    out_file = '%s/chr.%s.%d.famgen.counts.txt' % (args.out_dir, args.chrom, batch_num)
    print('saving to %s' % out_file)

    pos_data = np.load('%s/chr.%s.%d.gen.coordinates.npy' % (args.data_dir, args.chrom, batch_num))
    if pos_data.shape[0]>0:
        is_snp = pos_data[:, 2].astype(bool)
        is_pass = pos_data[:, 3].astype(bool)

        is_ok_include = np.ones(is_snp.shape, dtype=bool)
        if include_regions is not None:
            insert_loc = np.searchsorted(include_regions, pos_data[:, 1])
            is_ok_include = np.remainder(insert_loc, 2)==1

        is_ok_exclude = np.ones(is_snp.shape, dtype=bool)
        if exclude_regions is not None:
            insert_loc = np.searchsorted(exclude_regions, pos_data[:, 1])
            is_ok_exclude = np.remainder(insert_loc, 2)==0

        print('not SNP', np.sum(~is_snp))
        print('not PASS', np.sum(~is_pass))
        print('filtered by include', np.sum(~is_ok_include))
        print('filtered by exclude', np.sum(~is_ok_exclude))

        # Pull data together
        A = sparse.load_npz('%s/%s' % (args.data_dir, gen_file))

        # filter out snps
        A = A[:, is_snp & is_pass & is_ok_include & is_ok_exclude]
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
            
