import numpy as np
from collections import defaultdict, Counter, namedtuple
from itertools import product
import scipy.stats
import json
from os import listdir

BaselineCounts = namedtuple('BaselineCounts', ['counts', 'samples', 'families', 'family_sizes', 'is_child', 'is_mom', 'is_dad'])
Samples = namedtuple('Samples', ['sample_ids', 'families', 'family_sizes', 'is_child', 'is_mom', 'is_dad'])
PrecisionRecall = namedtuple('PrecisionRecall', ['precision1', 'recall1', 'F11', 'precision2', 'recall2', 'F12',
                    'precision1_upper_bound', 'recall1_upper_bound', 'F11_upper_bound', 'precision2_upper_bound', 'recall2_upper_bound', 'F12_upper_bound',
                    'error_rate', 'error_rate_no_missing'])

def pull_samples(data_dir, chroms, dot_in_name=False):
    sample_to_chroms = defaultdict(set)
    sample_to_family = dict()
    family_to_size = dict()
    children, moms, dads = set(), set(), set()
    chrom_batches = set()
    
    # pull counts from famgen
    for i, chrom in enumerate(chroms):
        print(chrom, end=' ')

        famgen_files = sorted([f for f in listdir(data_dir) if ('chr.%s.' % chrom) in f and '.famgen.counts.txt' in f], key=lambda x: int(x.split('.')[2]))
        for famgen_file in famgen_files:
            batch_num = int(famgen_file.split('.')[2])

            with open('%s/%s' % (data_dir, famgen_file), 'r') as f:
                for line in f:
                    pieces = line.strip().split('\t')
                    famkey, inds = pieces[:2]
                    #famkey = famkey.split('.')[0]

                    if dot_in_name:
                        # unfortunately, ssc uses . in their sample names
                        inds = inds.split('.')
                        inds = ['%s.%s' % (inds[i], inds[i+1]) for i in range(0, len(inds), 2)]
                    else:
                        inds = inds.split('.')
                        
                    m = len(inds)
                    family_to_size[famkey] = m
                    for ind in inds:
                        sample_to_chroms[ind].add((chrom, batch_num))
                        sample_to_family[ind] = famkey
                        chrom_batches.add((chrom, batch_num))

                    moms.add(inds[0])
                    dads.add(inds[1])
                    children.update(inds[2:])    
            
    multigen = children & (moms | dads)
    print('\nRemoving %d individuals involved in multiple generations' % len(multigen))
    
    children = children - multigen
    moms = moms - multigen
    dads = dads - multigen
    
    samples = sorted(children | moms | dads)
    families = [sample_to_family[x] for x in samples]
    family_sizes = np.array([family_to_size[x] for x in families])
    is_child = np.array([x in children for x in samples])
    is_mom = np.array([x in moms for x in samples])
    is_dad = np.array([x in dads for x in samples])
    
    return Samples(samples, families, family_sizes, is_child, is_mom, is_dad)

def pull_error_rates(samples, param_file, gens, obss):
    
    with open(param_file, 'r') as f:
        params = json.load(f)
    
    rates = np.zeros((len(samples.sample_ids), len(gens), len(obss)))
    rates[:] = np.nan
                
    for i, (sample_id, family) in enumerate(zip(samples.sample_ids, samples.families)):
        key = '%s.%s' % (family, sample_id)
        if key in params:
            rates[i, :, :] = [[params[key]['-log10(P[obs=%s|true_gen=%s])' % (o, g)] for o in obss] for g in gens]
        
    return rates

def pull_error_counts(samples, param_file, gens, obss):

    with open(param_file, 'r') as f:
        params = json.load(f)
    
    counts = np.zeros((len(samples.sample_ids), len(gens), len(obss)))
    counts[:] = np.nan
                
    for i, (sample_id, family) in enumerate(zip(samples.sample_ids, samples.families)):
        key = '%s.%s' % (family, sample_id)
        if key in params:
            counts[i, :, :] = [[params[key]['E[obs=%s, true_gen=%s]' % (o, g)] for o in obss] for g in gens]
        
    return counts

def pull_num_sites(d, chroms):
    all_num_sites = []
    for famgen_dir, param_file in zip(d['family_genotype_dir'], d['param_file']):
        samples = pull_samples(famgen_dir, chroms)

        with open(param_file, 'r') as f:
            params = json.load(f)
        
        num_sites = np.zeros((len(samples.sample_ids),))
        num_sites[:] = np.nan
                    
        for i, (sample_id, family) in enumerate(zip(samples.sample_ids, samples.families)):
            key = '%s.%s' % (family, sample_id)
            if key in params:
                num_sites[i] = sum([params[key]['observed_%s' % g] for g in ['0/0', '0/1', '1/1', './.']])
        all_num_sites.append(num_sites)
        
    return np.hstack(all_num_sites)

def pull_precision_recall(samples, param_file):
    with open(param_file, 'r') as f:
        params = json.load(f)

    het_precision = np.zeros((len(samples.sample_ids),))
    het_precision[:] = np.nan
    het_recall = np.zeros((len(samples.sample_ids),))
    het_recall[:] = np.nan
    het_F1 = np.zeros((len(samples.sample_ids),))
    het_F1[:] = np.nan

    homalt_precision = np.zeros((len(samples.sample_ids),))
    homalt_precision[:] = np.nan
    homalt_recall = np.zeros((len(samples.sample_ids),))
    homalt_recall[:] = np.nan
    homalt_F1 = np.zeros((len(samples.sample_ids),))
    homalt_F1[:] = np.nan

    het_precision_upper_bound = np.zeros((len(samples.sample_ids),))
    het_precision_upper_bound[:] = np.nan
    het_recall_upper_bound = np.zeros((len(samples.sample_ids),))
    het_recall_upper_bound[:] = np.nan
    het_F1_upper_bound = np.zeros((len(samples.sample_ids),))
    het_F1_upper_bound[:] = np.nan

    homalt_precision_upper_bound = np.zeros((len(samples.sample_ids),))
    homalt_precision_upper_bound[:] = np.nan
    homalt_recall_upper_bound = np.zeros((len(samples.sample_ids),))
    homalt_recall_upper_bound[:] = np.nan
    homalt_F1_upper_bound = np.zeros((len(samples.sample_ids),))
    homalt_F1_upper_bound[:] = np.nan

    error_rate = np.zeros((len(samples.sample_ids),))
    error_rate[:] = np.nan
    error_rate_no_missing = np.zeros((len(samples.sample_ids),))
    error_rate_no_missing[:] = np.nan
                
    for i, (sample_id, family) in enumerate(zip(samples.sample_ids, samples.families)):
        key = '%s.%s' % (family, sample_id)
        if key in params:
            het_precision[i] = params[key]['precision_0/1']
            het_recall[i] = params[key]['recall_0/1']
            het_F1[i] = params[key]['F1_0/1']
            homalt_precision[i] = params[key]['precision_1/1']
            homalt_recall[i] = params[key]['recall_1/1']
            homalt_F1[i] = params[key]['F1_1/1']

            het_precision_upper_bound[i] = params[key]['upper_bound[precision_0/1]']
            het_recall_upper_bound[i] = params[key]['upper_bound[recall_0/1]']
            het_F1_upper_bound[i] = params[key]['upper_bound[F1_0/1]']
            homalt_precision_upper_bound[i] = params[key]['upper_bound[precision_1/1]']
            homalt_recall_upper_bound[i] = params[key]['upper_bound[recall_1/1]']
            homalt_F1_upper_bound[i] = params[key]['upper_bound[F1_1/1]']

            match = params[key]['E[obs=0/0, true_gen=0/0]'] + params[key]['E[obs=0/1, true_gen=0/1]'] + params[key]['E[obs=1/1, true_gen=1/1]'] 
            error = params[key]['E[obs=0/1, true_gen=0/0]'] + params[key]['E[obs=1/1, true_gen=0/0]'] + params[key]['E[obs=./., true_gen=0/0]'] + \
            		params[key]['E[obs=0/0, true_gen=0/1]'] + params[key]['E[obs=1/1, true_gen=0/1]'] + params[key]['E[obs=./., true_gen=0/1]'] + \
            		params[key]['E[obs=0/0, true_gen=1/1]'] + params[key]['E[obs=0/1, true_gen=1/1]'] + params[key]['E[obs=./., true_gen=1/1]']
            error_no_missing = params[key]['E[obs=0/1, true_gen=0/0]'] + params[key]['E[obs=1/1, true_gen=0/0]'] + \
            		params[key]['E[obs=0/0, true_gen=0/1]'] + params[key]['E[obs=1/1, true_gen=0/1]'] + \
            		params[key]['E[obs=0/0, true_gen=1/1]'] + params[key]['E[obs=0/1, true_gen=1/1]']

            error_rate[i] = error/(error+match)
            error_rate_no_missing[i] = error_no_missing/(error_no_missing+match)
        
    return PrecisionRecall(het_precision, het_recall, het_F1, homalt_precision, homalt_recall, homalt_F1,
        het_precision_upper_bound, het_recall_upper_bound, het_F1_upper_bound, homalt_precision_upper_bound, homalt_recall_upper_bound, homalt_F1_upper_bound,
        error_rate, error_rate_no_missing)

def mad_outlier_detection(x):
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    return ~np.isnan(x) & (0.6745*(x - med)/mad > 3.5)

def process_datasets(datasets, chroms):
    for d in datasets:
        print(d)
        if not isinstance(d['family_genotype_dir'], list):
            d['family_genotype_dir'] = [d['family_genotype_dir']]
            d['param_file'] = [d['param_file']]
            
        for key in ['precision1', 'precision2', 'recall1', 'recall2', 'F11', 'F12']:
            d['%s_estimates' % key] = []
            d['%s_upper_bounds' % key] = []
            d['%s_samples' % key] = []
        d['error_rate'] = []
        d['error_rate_no_missing'] = []
        
        for famgen_dir, param_file in zip(d['family_genotype_dir'], d['param_file']):
        
            samples = pull_samples(famgen_dir, chroms)
            include_sample = d['sample_filter'](samples)
            precision_recall = pull_precision_recall(samples, param_file) 
            
            # pull num missing
            with open(param_file, 'r') as f:
                params = json.load(f)

            family_size = np.zeros((len(samples.sample_ids),))
            family_size[:] = np.nan
            nm = np.zeros((len(samples.sample_ids),))
            nm[:] = np.nan
            missing = np.zeros((len(samples.sample_ids),))
            missing[:] = np.nan

            for i, (sample_id, family, include) in enumerate(zip(samples.sample_ids, samples.families, include_sample)):
                key = '%s.%s' % (family, sample_id)
                if include and (key in params):
                    family_size[i] = params[key]['family_size']
                    nm[i] = np.log10(params[key]['nonmendelian_count'])
                    missing[i] = np.log10(params[key]['missing_count'])
            
            is_nm_outlier = mad_outlier_detection(nm)
            is_missing_outlier = mad_outlier_detection(missing)
            print('num removed due to filter', np.sum(~include_sample))
            print('num removed due to missing', np.sum(is_missing_outlier))
            print('num removed due to too many non-Mendelian sites', np.sum(is_nm_outlier))


            for key in ['precision1', 'precision2', 'recall1', 'recall2', 'F11', 'F12']:
                estimates = getattr(precision_recall, key)
                upper_bounds = getattr(precision_recall, '%s_upper_bound' % key)
                indices = include_sample & ~np.isnan(estimates) & ~np.isnan(upper_bounds) & ~is_nm_outlier & ~is_missing_outlier
                d['%s_estimates' % key].append(estimates[indices])
                d['%s_upper_bounds' % key].append(upper_bounds[indices])
                d['%s_samples' % key].extend([samples.sample_ids[i] for i in np.where(indices)[0]])
            d['error_rate'].extend(precision_recall.error_rate[indices])
            d['error_rate_no_missing'].extend(precision_recall.error_rate_no_missing[indices])
            print(d['precision1_estimates'][-1].shape)
                
        for key in ['precision1', 'precision2', 'recall1', 'recall2', 'F11', 'F12']:
            d['%s_estimates' % key] = np.hstack(d['%s_estimates' % key])
            d['%s_upper_bounds' % key] = np.hstack(d['%s_upper_bounds' % key])
            print(d['%s_estimates' % key].shape, d['%s_upper_bounds' % key].shape)
        d['error_rate'] = np.hstack(d['error_rate'])
        d['error_rate_no_missing'] = np.hstack(d['error_rate_no_missing'])



