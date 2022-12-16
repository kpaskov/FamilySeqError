import numpy as np

## ------------------------------------ Calculate Various Metrics ------------------------------------

def add_observed_counts(params, counts, j, m, gens, obss):
    for i, obs in enumerate(obss):
        params['observed_%s' % obs] = int(np.sum(counts[tuple(i if x==j else slice(None, None, None) for x in range(m))]))

def add_estimated_error_rates(params, error_rates, lower_bounds, j, gens, obss):
    for gen_index, gen in enumerate(gens):
        for obs_index, obs in enumerate(obss):
            params['-log10(P[obs=%s|true_gen=%s])' % (obs, gen)] = float(-np.log10(error_rates[j, gen_index, obs_index]))
            params['lower_bound[-log10(P[obs=%s|true_gen=%s])]' % (obs, gen)] = float(-np.log10(lower_bounds[j, gen_index, obs_index]))

def add_expected_counts(params, gens, obss):
    # we assume error rates are low, so the number of times we observe a genotype is a good estimate of the number of times this genotype actually occurs.
    for gen_index, gen in enumerate(gens):
        for obs_index, obs in enumerate(obss):
            params['E[obs=%s, true_gen=%s]' % (obs, gen)] = params['observed_%s' % gen] * (10.0**-params['-log10(P[obs=%s|true_gen=%s])' % (obs, gen)])
            params['lower_bound[E[obs=%s, true_gen=%s]]' % (obs, gen)] = params['observed_%s' % gen] * (10.0**-params['lower_bound[-log10(P[obs=%s|true_gen=%s])]' % (obs, gen)])

def add_precision_recall(params, gens, obss):
    # precision: TP/(TP + FP)
    # let n_0 = # of times the real genotype is 0/0
    # E[TP] = n_1 * p_11
    # E[FP] = n_0 * p_01

    # we again assume error rates are low, so the number of times we observe a genotype is a good estimate of the number of times this genotype actually occurs.

    for var in gens:
        TP = params['E[obs=%s, true_gen=%s]' % (var, var)]
        FP = np.sum([params['E[obs=%s, true_gen=%s]' % (var, gen)] for gen in gens if var != gen])
        FN = np.sum([params['E[obs=%s, true_gen=%s]' % (obs, var)] for obs in obss if var != obs])
        

        FP_lb = np.sum([params['lower_bound[E[obs=%s, true_gen=%s]]' % (var, gen)] for gen in gens if var != gen])
        FN_lb = np.sum([params['lower_bound[E[obs=%s, true_gen=%s]]' % (obs, var)] for obs in obss if var != obs])

        params['precision_%s' % var] = TP/(TP+FP)
        params['recall_%s' % var] = TP/(TP+FN)
        params['FPR_%s' % var] = FP/(TP+FP)
        params['FNR_%s' % var] = FN/(TP+FN)
        params['F1_%s' % var] = TP/(TP + 0.5*(FP+FN))

        params['upper_bound[precision_%s]' % var] = TP/(TP+FP_lb)
        params['upper_bound[recall_%s]' % var] = TP/(TP+FN_lb)
        params['upper_bound[F1_%s]' % var] = TP/(TP + 0.5*(FP_lb+FN_lb))