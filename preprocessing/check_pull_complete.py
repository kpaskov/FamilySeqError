import argparse
from pysam import VariantFile
import os.path
import math

parser = argparse.ArgumentParser(description='Pull genotypes.')
parser.add_argument('vcf_file', type=str, help='VCF file to pull from.')
parser.add_argument('data_dir', type=str, help='Directory containing genotype data.')
parser.add_argument('file_type', type=str, help='Type of file to check.')
parser.add_argument('--batch_size', type=int, default=None, help='Restrict number of positions per file to batch_size.')
args = parser.parse_args()

chroms = [str(x) for x in range(1, 23)]

def process_header(vcf):
    return sample_ids, vcf.header.contigs

vcf = VariantFile(args.vcf_file)
contigs = vcf.header.contigs

missing_files = []
for chrom in chroms:
	# pull contig length
	contig = None
	if chrom in contigs:
	    contig = contigs[chrom]
	elif 'chr%s' % chrom in contigs:
	    contig = contigs['chr%s' % chrom]
	else:
	    raise Exception('Trouble finding contig', chrom, 'in', contig_names)
	print('Chrom %s length %d' % (chrom, contig.length))


	if args.batch_size is None:
		num_batches = 1
	else:
		num_batches = math.ceil(contig.length/args.batch_size)
		

	for batch_num in range(num_batches):
		filename = '%s/chr.%s.%d.%s' % (args.data_dir, chrom, batch_num, args.file_type)
		if not os.path.isfile(filename):
			missing_files.append(filename)
print('Missing files:', missing_files)


