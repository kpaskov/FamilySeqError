import argparse
import os.path
import math
import json

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('data_dir', type=str, help='Directory containing genotype data.')
parser.add_argument('file_type', type=str, help='Type of file to check.')
args = parser.parse_args()

chroms = [str(x) for x in range(1, 23)]

# grab assembly
with open('%s/info.json' % args.data_dir, 'r') as f:
	info = json.load(f)
	assembly = info['assembly']
	batch_size = info['batch_size']

with open('data/chrom_lengths%s.json' % assembly, 'r') as f:
	chrom_lengths = json.load(f)

missing_files = []
for chrom in chroms:
	chrom_length = chrom_lengths[chrom]

	if batch_size is None:
		num_batches = 1
	else:
		num_batches = math.ceil(chrom_length/batch_size)
		

	for batch_num in range(num_batches):
		filename = '%s/chr.%s.%d.%s' % (args.data_dir, chrom, batch_num, args.file_type)
		if not os.path.isfile(filename):
			missing_files.append(filename)
print('Missing files:', missing_files)


