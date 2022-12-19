import sys
import numpy as np
import scipy.stats
from itertools import product
import cvxpy as cp
from collections import Counter, defaultdict
import json
import os
import random
from calculate_metrics import add_observed_counts, add_estimated_error_rates, add_expected_counts, add_precision_recall

import argparse

parser = argparse.ArgumentParser(description='Estimate parameters.')
parser.add_argument('data_dir', type=str, help='Family genotype count directory.')
parser.add_argument('ped_file', type=str, help='Ped file.')
parser.add_argument('--is_ngs', action='store_true', default=False, help='True if this data is NGS. The important point is whether or not sites where all individuals are homozygous reference are sometimes dropped from the VCF. If this happens, use flag --is_ngs')
parser.add_argument('--family', type=str, default=None, help='Estimate parameters for a given family only.')
parser.add_argument('--count_type', type=str, default=None, help='Name of count type. Used to differentiate between counts in high-complexity vs low-complexity regions, for example.')
args = parser.parse_args()

if not os.path.exists('%s/sequencing_error_rates' % args.data_dir):
    os.makedirs('%s/sequencing_error_rates' % args.data_dir)

if args.count_type is None:
    input_dir = '%s/family_genotype_counts' % args.data_dir
    out_file = '%s/sequencing_error_rates/errors.json' % args.data_dir
else:
    input_dir = '%s/family_genotype_counts/%s' % (args.data_dir, args.count_type)
    out_file = '%s/sequencing_error_rates/%s_errors.json' % (args.data_dir, args.count_type)

chroms = ['X']

gens = ['0/0', '0/1', '1/1'] # homref, het, homalt

obss = ['0/0', '0/1', '1/1', './.']

errors = [
    ('0/0', '0/1'),
    ('0/0', '1/1'), 
    ('0/1', '0/0'), 
    ('0/1', '1/1'), 
    ('1/1', '0/0'), 
    ('1/1', '0/1'), 
    ('0/0', './.'),
    ('0/1', './.'),
    ('1/1', './.'),
]

father_errors = [
    ('0/0', '0/1'),
    ('0/0', '1/1'),
    ('1/1', '0/0'),
    ('1/1', '0/1'),
    ('0/0', './.'),
    ('1/1', './.'),
]
father_error_to_index = dict([(x, errors.index(x)) for x in father_errors])

mother_errors = [
    ('0/0', '0/1'),
    ('0/0', '1/1'), 
    ('0/1', '0/0'), 
    ('0/1', '1/1'), 
    ('1/1', '0/0'), 
    ('1/1', '0/1'), 
    ('0/0', './.'),
    ('0/1', './.'),
    ('1/1', './.'),
]
mother_error_to_index = dict([(x, errors.index(x)) for x in mother_errors])

son_errors = [
    ('0/0', '0/1'),
    ('0/0', '1/1'),
    ('1/1', '0/0'),
    ('1/1', '0/1'),
    ('0/0', './.'),
    ('1/1', './.'),
]
son_error_to_index = dict([(x, errors.index(x)) for x in son_errors])

daughter_errors = [
    ('0/0', '0/1'),
    ('0/0', '1/1'), 
    ('0/1', '0/0'), 
    ('0/1', '1/1'), 
    ('1/1', '0/0'), 
    ('1/1', '0/1'), 
    ('0/0', './.'),
    ('0/1', './.'),
    ('1/1', './.'),
]
daughter_error_to_index = dict([(x, errors.index(x)) for x in daughter_errors])

son_mendelian_trios = {
    ('0/0', '0/0', '0/0'),
    ('0/0', '1/1', '0/0'),
    ('0/1', '0/0', '0/0'), ('0/1', '0/0', '1/1'),
    ('0/1', '1/1', '0/0'), ('0/1', '1/1', '1/1'),
    ('1/1', '0/0', '1/1'),
    ('1/1', '1/1', '1/1')
}

daughter_mendelian_trios = {
    ('0/0', '0/0', '0/0'),
    ('0/0', '1/1', '0/1'),
    ('0/1', '0/0', '0/0'), ('0/1', '0/0', '0/1'),
    ('0/1', '1/1', '0/1'), ('0/1', '1/1', '1/1'),
    ('1/1', '0/0', '0/1'),
    ('1/1', '1/1', '1/1')
}

print('num error types', len(errors))



# ------------------------------------ Pull Data ------------------------------------

# Sex (1=male; 2=female; other=unknown)
sample_id_to_sex = dict()
with open(args.ped_file, 'r') as f:
    for line in f:
        pieces = line.strip().split('\t')
        if len(pieces) >= 6:
            fam_id, sample_id, f_id, m_id, sex, disease_status = pieces[0:6]
            sample_id_to_sex[sample_id] = sex

family_to_counts = dict()
family_to_inds = dict()
family_to_indices = dict()
for i, chrom in enumerate(chroms):
    print(chrom, end=' ')

    count_files = sorted([f for f in os.listdir(args.data_dir) if ('chr.%s.' % chrom) in f and 'famgen.counts.txt' in f])
    for count_file in count_files:
        with open('%s/%s' % (args.data_dir, count_file), 'r') as f:
            for line in f:
                famkey = line.strip().split('\t', maxsplit=1)[0]
                if args.family==None or args.family==famkey:
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


famkeys = sorted(family_to_inds.keys())
print('Families', len(famkeys))
print('Families of each size', Counter([len(inds) for fkey, inds in family_to_inds.items()]))



# ------------------------------------ Poisson Regression ------------------------------------

def get_mendelian(m, inds):
    is_male = [sample_id_to_sex[x]=='1' for x in inds]

    # differentiate mendelian and non-mendelian famgens
    is_mendelian = np.ones((len(obss),)*m, dtype=bool)
    for famgen in product(range(len(obss)), repeat=m):
        is_mend = True
        for j in range(2, m):
            if is_male[j]:
                if not (tuple(obss[famgen[x]] for x in [0, 1, j]) in son_mendelian_trios):
                    is_mend = False
            else:
                if not (tuple(obss[famgen[x]] for x in [0, 1, j]) in daughter_mendelian_trios):
                    is_mend = False
        is_mendelian[famgen] = is_mend
    return is_mendelian

def has_variant(x):
    return len([y for y in x if y>0])>0

def estimate_error_rates(is_mendelian, error_to_index, obs_counts, true_counts):

    # -------------------- set up problem --------------------
    m = len(obs_counts.shape)
    #famgens = [x for x in product(np.arange(len(obss)), repeat=m)]
    famgens = [x for x in zip(*np.where(~is_mendelian))]

    if args.is_ngs:
        # if we're working with NGS data, we don't know the real counts of famgens without variants
        # because they may have been excluded from the vcf
        famgens = [x for x in famgens if has_variant(x)]

    m = len(obs_counts.shape)
    X = np.zeros((len(famgens), len(errors)*m), dtype=float)
    y = np.zeros((len(famgens),), dtype=int)

    for k, fg in enumerate(famgens):
        for i, j in product(range(4), range(m)):
            error = (obss[i], obss[fg[j]])
            if error in error_to_index[j]:
                neighbor = tuple(i if k==j else fg[k] for k in range(m))
                if is_mendelian[neighbor]:
                    error_index = error_to_index[j][error] + j*len(errors)
                    X[k, error_index] += true_counts[neighbor]
        y[k] = obs_counts[fg]
    X[np.any(np.isnan(X), axis=1), :] = 0

    #print([famgens[i] for i in np.where(X[:, error_to_index[('0/0', '0/1')]+2*len(errors)] != 0)[0]])


    print('Removing zero rows:', np.sum(np.sum(X, axis=1)==0))
    indices = np.where((np.sum(X, axis=1) != 0))[0]
    X = X[indices, :]
    y = y[indices]

    is_zero = np.sum(X, axis=0)==0
    print('Removing zero cols:', [(np.floor(i/len(errors)), errors[i % len(errors)]) for i in np.where(is_zero)[0]])
    X = X[:, ~is_zero]
    old_col_index_to_new = dict([(old_index, new_index) for new_index, old_index in enumerate(np.where(~is_zero)[0])])



    print(X.shape, y.shape)

    # -------------------- solve problem --------------------

    print('Estimating...', X.shape, y.shape)
    alpha = 1.0/np.max(X)

    # cvxpy
    n = cp.Variable(X.shape[1])
    mu = np.sum(X, axis=0)
    objective = cp.Minimize(alpha*mu@n - alpha*y@cp.log(X@n))

    # Wilson score interval so that if we don't observe any errors, then we take the 95% confidence interval
    z = 1.96
    lower_bound = ((z*z)/2)/(mu+(z*z))
    constraints = [n >= 0, n <= 1]
    prob = cp.Problem(objective, constraints)

    result = prob.solve(solver='ECOS', max_iters=10000, verbose=False)
    print(prob.status)

    #print(n.value, n.value.shape)
    ns = np.asarray([v for v in n.value])

    if prob.status != 'optimal' and prob.status != 'optimal_inaccurate':
        raise Error('Parameters not fully estimated.')

    # -------------------- reformat solution --------------------

    error_rates = np.zeros((len(inds), len(gens), len(obss)), dtype=float)
    lower_bounds = np.zeros((len(inds), len(gens), len(obss)), dtype=float)
    error_rates[:] = np.nan
    lower_bounds[:] = np.nan
    for k in range(len(errors)*len(inds)):
        if k in old_col_index_to_new:
            e = errors[k%len(errors)]
            ind_index = int(np.floor(k/len(errors)))
            error_rates[ind_index, gens.index(e[0]), obss.index(e[1])] = ns[old_col_index_to_new[k]]
            lower_bounds[ind_index, gens.index(e[0]), obss.index(e[1])] = lower_bound[old_col_index_to_new[k]]

    # now fill in P(obs=true_gen)
    for i, gen in enumerate(gens):
        if gen in obss:
            error_rates[:, i, obss.index(gen)] = 1-np.sum(error_rates[:, i, [k for k in range(len(obss)) if k != i]], axis=1)

    return error_rates, lower_bounds, np.sum([obs_counts[fg] for fg in famgens if has_variant(fg)])


# ------------------------------------ Estimate Error Rates ------------------------------------

params = {}
baseline_match = {(0, 0), (1, 1), (2, 2)}

num_error_families = 0
for i, famkey in enumerate(famkeys):
    print(famkey)
    #try:
    inds = family_to_inds[famkey]
    m = len(inds)

    is_mendelian = get_mendelian(m, inds)
    error_to_index = [mother_error_to_index, father_error_to_index] + \
    [son_error_to_index if sample_id_to_sex[x]=='1' else daughter_error_to_index for x in inds[2:]]
    #error_to_index = [parent_error_to_index]*2 + [child_error_to_index]*(m-2)
    #error_to_index = [set()]*2 + [child_error_to_index]*(m-2)
    #error_to_index = [child_error_to_index]*m

    counts = family_to_counts[famkey]
    real_counts = counts.copy().astype(float)


    print('-------------Estimate error rates-------------')
    error_rates, lower_bounds, nm_count = estimate_error_rates(is_mendelian, error_to_index, counts, real_counts)

    for j in range(len(inds)):
        # observed counts
        ind_params = {'nonmendelian_count': int(nm_count), 
                        'family_size': int(m), 
                        'missing_count': int(np.sum(counts.take(indices=3, axis=j))),
                        'total_count': int(np.sum(counts))}
        add_observed_counts(ind_params, counts, j, m, gens, obss)
        add_estimated_error_rates(ind_params, error_rates, lower_bounds, j, gens, obss)
        add_expected_counts(ind_params, gens, obss)
        add_precision_recall(ind_params, gens, obss)

        params[famkey[0] + '.' + inds[j]] = ind_params

    #except Exception as err:
    #    num_error_families += 1
    #    print('ERROR', err)

print('Total errors', num_error_families)

# ------------------------------------ Write to file ------------------------------------

with open(out_file, 'w+') as f:
    json.dump(params, f, indent=4)
