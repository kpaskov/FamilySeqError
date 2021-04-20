import numpy as np
from collections import defaultdict, Counter, namedtuple
from itertools import product
import scipy.stats
import json
from os import listdir
from compare_rates import pull_samples, pull_error_rates
from calculate_metrics import add_observed_counts, add_estimated_error_rates, add_expected_counts, add_precision_recall

import argparse

parser = argparse.ArgumentParser(description='Estimate parameters.')
parser.add_argument('param_file', type=str, help='Parameter file for families you want to use to find the median.')
parser.add_argument('data_dir', type=str, help='Family genotype count directory for families you want to use to find the median.')
parser.add_argument('new_data_dir', type=str, help='Family genotype count directory for families you want to create new param file for.')
parser.add_argument('out_file', type=str, help='Output file.')
args = parser.parse_args()

chroms = [str(x) for x in range(1, 23)]
gens = ['0/0', '0/1', '1/1']
obss = ['0/0', '0/1', '1/1', './.']

# get median error rates
samples = pull_samples(args.data_dir, [str(x) for x in range(1, 23)])
rates = pull_error_rates(samples, args.param_file, gens, obss)

child_rates = rates[samples.is_child, :, :]
print('nans', np.sum(np.isnan(child_rates)))

median_rates = 10**-np.nanmedian(child_rates, axis=0)

# readjust rates to sum to 1
for i, gen in enumerate(gens):
	true_obs_index = obss.index(gen)
	error_rate = np.sum(median_rates[i, :]) - median_rates[i, true_obs_index]
	median_rates[i, true_obs_index] = 1-error_rate

assert np.all(median_rates>0) and np.all(median_rates<1)
print('median error rates')
print(median_rates)


# ------------------------------------ Pull Data ------------------------------------

family_to_counts = dict()
family_to_inds = dict()
family_to_indices = dict()
for i, chrom in enumerate(chroms):
    print(chrom, end=' ')
    
    count_files = sorted([f for f in listdir(args.new_data_dir) if ('chr.%s.' % chrom) in f and 'famgen.counts.txt' in f])
    for count_file in count_files:
        with open('%s/%s' % (args.new_data_dir, count_file), 'r') as f:
            for line in f:
                famkey = line.strip().split('\t', maxsplit=1)[0]
                pieces = line.strip().split('\t')
                inds = pieces[1].split('.')
                m = len(inds)

                counts = np.zeros((len(obss),)*m, dtype=int)
                for g, c in zip(product(range(len(obss)), repeat=m), pieces[2:]):
                    counts[g] = int(c)

                if (famkey, inds[0], inds[1]) in family_to_inds:
                    assert family_to_inds[(famkey, inds[0], inds[1])] == inds
                    counts += family_to_counts[(famkey, inds[0], inds[1])]
                else:
                    family_to_inds[(famkey, inds[0], inds[1])] = inds
                    family_to_indices[(famkey, inds[0], inds[1])] = list(range(m))
                    
                family_to_counts[(famkey, inds[0], inds[1])] = counts


famkeys = [x for x in sorted(family_to_inds.keys()) if len(family_to_inds[x])==3]
print('Families', len(famkeys))
print('Families of each size', Counter([len(inds) for fkey, inds in family_to_inds.items()]))

params = {}
for i, famkey in enumerate(famkeys):
    print(famkey)
    inds = family_to_inds[famkey]
    m = len(inds)
    
    counts = family_to_counts[famkey]
        
    print('-------------Estimate error rates-------------')
    error_rates = np.vstack([median_rates[np.newaxis, :, :] for _ in range(m)])
    print(error_rates.shape)
    lower_bounds = np.zeros((m, len(gens), len(obss)), dtype=float)
    lower_bounds[:] = np.nan
    
    for j in range(len(inds)):
        # observed counts
        ind_params = {'nonmendelian_count': -1, 
                        'family_size': int(m), 
                        'missing_count': int(np.sum(counts.take(indices=3, axis=j))),
                        'total_count': int(np.sum(counts))}
        add_observed_counts(ind_params, counts, j, m, gens, obss)
        add_estimated_error_rates(ind_params, error_rates, lower_bounds, j, gens, obss)
        add_expected_counts(ind_params, gens, obss)
        add_precision_recall(ind_params, gens, obss)

        params[famkey[0] + '.' + inds[j]] = ind_params

# ------------------------------------ Write to file ------------------------------------

with open(args.out_file, 'w+') as f:
    json.dump(params, f, indent=4)





