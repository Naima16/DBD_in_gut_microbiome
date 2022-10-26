import numpy as np
import bz2
from collections import defaultdict
import sys

species = sys.argv[1]
genes_dir = '/u/home/d/daisyche/dbd/data/genes'
fname = '%s/%s/genes_presabs.txt.bz2' % (genes_dir, species)
file = bz2.BZ2File(fname, 'r')

# Total number of present genes in each sample
num_pres_genes = defaultdict(int)
total_genes = 0

# Total number of present 99% core genes in each sample
num_pres_core99_genes = defaultdict(int)
total_core99_genes = 0

# Total number of present 95% core genes in each sample
num_pres_core95_genes = defaultdict(int)
total_core95_genes = 0

# Total number of present 90% core genes in each sample
num_pres_core90_genes = defaultdict(int)
total_core90_genes = 0

# Samples
samples = file.readline().strip().split('\t')[1:]
num_samples = len(samples)

# Loop over genes
for line in file:
	is_pres_arr = map(lambda x: int(x), line.strip().split('\t')[1:])
	gene_prev = sum(is_pres_arr)
	
	total_genes += 1
	if (gene_prev/float(num_samples)) > 0.99:
		total_core99_genes += 1
	if (gene_prev/float(num_samples)) > 0.95:
		total_core95_genes += 1
	if (gene_prev/float(num_samples)) > 0.9:
		total_core90_genes += 1
	
	for i in range(len(is_pres_arr)):
		sample, is_pres = samples[i], is_pres_arr[i]
		
		num_pres_genes[sample] += is_pres
		if (gene_prev/float(num_samples)) > 0.99:
			num_pres_core99_genes[sample] += is_pres
		if (gene_prev/float(num_samples)) > 0.95:
			num_pres_core95_genes[sample] += is_pres
		if (gene_prev/float(num_samples)) > 0.9:
			num_pres_core90_genes[sample] += is_pres

# Accessory genes = total # genes - (# core genes, which are present in xx% of samples)

num_pres_acc99_genes = {sample: (num_pres_genes[sample] - num_pres_core99_genes[sample]) for sample in samples}
total_acc99_genes = total_genes - total_core99_genes

num_pres_acc95_genes = {sample: (num_pres_genes[sample] - num_pres_core95_genes[sample]) for sample in samples}
total_acc95_genes = total_genes - total_core95_genes

num_pres_acc90_genes = {sample: (num_pres_genes[sample] - num_pres_core90_genes[sample]) for sample in samples}
total_acc90_genes = total_genes - total_core90_genes

# Average the number of present genes for all samples that have a particular species

pres_genes_arr = num_pres_genes.values()
pres_genes_arr_no0 = filter(lambda num: num > 0, pres_genes_arr)
avg_pres_genes = np.array(pres_genes_arr_no0).mean()

genes_counts_dir = '/u/home/d/daisyche/dbd/data/gene_counts'
output_fname = '%s/gene_counts_%s.csv' % (genes_counts_dir, species)
output_file = open(output_fname, 'w')
output_file.write(','.join(['species_name','sample_id','pres_genes','total_genes','pres_acc95_genes','total_acc95_genes','pres_acc90_genes','total_acc90_genes','species_avg_pres_genes']) + '\n')
for sample in samples:
	items = [species, sample, num_pres_genes[sample], total_genes, num_pres_acc95_genes[sample], total_acc95_genes, num_pres_acc90_genes[sample], total_acc90_genes, round(avg_pres_genes, 1)]
	output_file.write(','.join(map(lambda x: str(x), items)) + '\n')
