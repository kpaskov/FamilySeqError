import sys


data_dir = sys.argv[1]
chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']

for chrom in chroms:
    print(chrom, end=' ')

    with open('%s/chr.%s.gen.samples.txt' % (data_dir, chrom), 'r') as f:
        sample_ids = [x.strip() for x in f]

    with open('%s/chr.%s.gen.samples.txt' % (data_dir, chrom), 'w+') as f:
    	f.write('\n'.join([x.split('_')[-1] for x in sample_ids]))
