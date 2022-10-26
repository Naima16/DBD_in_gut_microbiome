from utils import config, parse_midas_data as pmd
import numpy as np
import bz2
from collections import defaultdict
import sys
import pickle

# Genes data directory
genes_dir = '/u/home/d/daisyche/dbd/data/genes'

# Output file
genes_counts_dir = '/u/home/d/daisyche/dbd/data/gene_counts'
output_fname = '%s/gene_counts_combined.csv' % (genes_counts_dir)
output_file = open(output_fname, 'w')
output_file.write(','.join(['species_name','sample_id','pres_genes','total_genes','community_pres_genes','community_pres_genes_nonredundant','pres_acc95_genes','total_acc95_genes','community_pres_acc95_genes','community_pres_acc95_genes_nonredundant','pres_acc90_genes','total_acc90_genes','community_pres_acc90_genes','community_pres_acc90_genes_nonredundant']) + '\n')

# Species list
species_list = pmd.parse_species_list()

# types we are considering
types = ['all', 'core99', 'core95', 'core90', 'acc99', 'acc95', 'acc90']

# type -> species -> sample -> set of present genes
pres_genes_dict = {type: {species: defaultdict(set) for species in species_list} for type in types}

# type -> species -> sample -> set of present gene clusters
pres_gene_clusters_dict = {type: {species: defaultdict(set) for species in species_list} for type in types}

# type -> sample -> set of present genes, aggregated across species
sample_pres_genes_dict = {type: defaultdict(set) for type in types}

# type -> species -> set of present gene clusters, aggregated across species
sample_pres_gene_clusters_dict = {type: defaultdict(set) for type in types}

# type -> species -> total genes considered (set)
total_genes_dict = {type: defaultdict(set) for type in types}

# get pickled gene cluster map
pdir = '%s/pickles' % config.data_directory
gene_cluster_map = pickle.load(open('%s/gene_cluster_map.pkl' % pdir, 'rb'))

for species in species_list:
	
	print("Working on %s..." % species)
	
	# Gene presence/absence info file
	fname = '%s/%s/genes_presabs.txt.bz2' % (genes_dir, species)
	file = bz2.BZ2File(fname, 'r')
	
	# Samples
	samples = file.readline().strip().split('\t')[1:]
	num_samples = len(samples)
	
	# Loop over genes
	for line in file:
		
		gene = line.strip().split('\t')[0]
		gene_cluster = gene_cluster_map[gene] if gene in gene_cluster_map else gene
		is_pres_arr = map(lambda x: int(x), line.strip().split('\t')[1:])
		gene_prev = sum(is_pres_arr)
		
		# Record in total_genes_dict
		total_genes_dict['all'][species].add(gene)
		
		if (gene_prev/float(num_samples)) > 0.99:
			total_genes_dict['core99'][species].add(gene)
		else:
			total_genes_dict['acc99'][species].add(gene)
		
		if (gene_prev/float(num_samples)) > 0.95:
			total_genes_dict['core95'][species].add(gene)
		else:
			total_genes_dict['acc95'][species].add(gene)
		
		if (gene_prev/float(num_samples)) > 0.9:
			total_genes_dict['core90'][species].add(gene)
		else:
			total_genes_dict['acc90'][species].add(gene)
		
		# Loop over samples
		for i in range(len(is_pres_arr)):
			sample, is_pres = samples[i], is_pres_arr[i] # is_pres: 0/1
			
			if is_pres == 1: # This gene is present in the sample
				
				# This is any gene in the species
				pres_genes_dict['all'][species][sample].add(gene)
				pres_gene_clusters_dict['all'][species][sample].add(gene_cluster)
				sample_pres_genes_dict['all'][sample].add(gene)
				sample_pres_gene_clusters_dict['all'][sample].add(gene_cluster)
				
				if (gene_prev/float(num_samples)) > 0.99: # This is a core95 gene for this species
					pres_genes_dict['core99'][species][sample].add(gene)
					pres_gene_clusters_dict['core99'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['core99'][sample].add(gene)
					sample_pres_gene_clusters_dict['core99'][sample].add(gene_cluster)
				else:
					pres_genes_dict['acc99'][species][sample].add(gene)
					pres_gene_clusters_dict['acc99'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['acc99'][sample].add(gene)
					sample_pres_gene_clusters_dict['acc99'][sample].add(gene_cluster)
				
				if (gene_prev/float(num_samples)) > 0.95: # This is a core95 gene for this species
					pres_genes_dict['core95'][species][sample].add(gene)
					pres_gene_clusters_dict['core95'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['core95'][sample].add(gene)
					sample_pres_gene_clusters_dict['core95'][sample].add(gene_cluster)
				else:
					pres_genes_dict['acc95'][species][sample].add(gene)
					pres_gene_clusters_dict['acc95'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['acc95'][sample].add(gene)
					sample_pres_gene_clusters_dict['acc95'][sample].add(gene_cluster)
				
				if (gene_prev/float(num_samples)) > 0.90: # This is a core90 gene for this species
					pres_genes_dict['core90'][species][sample].add(gene)
					pres_gene_clusters_dict['core90'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['core90'][sample].add(gene)
					sample_pres_gene_clusters_dict['core90'][sample].add(gene_cluster)
				else:
					pres_genes_dict['acc90'][species][sample].add(gene)
					pres_gene_clusters_dict['acc90'][species][sample].add(gene_cluster)
					sample_pres_genes_dict['acc90'][sample].add(gene)
					sample_pres_gene_clusters_dict['acc90'][sample].add(gene_cluster)

# Have a way to take set of genes, removing shared genes
from utils import midas_db_utils
shared_gene_dict = midas_db_utils.parse_midas_shared_genes_map()

# Designate gene cluster IDs (pick an arbitrary but consistent representative)
def get_gene_cluster_map():
	gene_cluster_map = {}
	covered_genes = []
	
	for gene in shared_gene_dict:
		if gene not in covered_genes: # not yet encountered in any previous cluster
			representative_gene = gene # designate key as representative gene
			
			gene_cluster_map[representative_gene] = representative_gene
			covered_genes.append(representative_gene)
			
			for shared_gene in shared_gene_dict[gene]:
				gene_cluster_map[shared_gene] = representative_gene
				covered_genes.append(shared_gene)
		
		# If encountered in previous cluster, it is already accounted for
	
	return gene_cluster_map

# This was too slow...
def get_nonredundant_genes(gene_set):
	nonredundant_genes = []
	redundant_genes = []
	
	for gene in gene_set:
		shared_genes = shared_gene_dict[gene]
		
		# Firstly, skip if this gene is redundant
		# with some previously traversed gene
		if gene in redundant_genes:
			continue
		
		# If this gene was not previously redundant,
		# assume it is nonredundant and
		# check if its shared genes are present
		nonredundant_genes.append(gene)
		
		for shared_gene in shared_genes:
			if shared_gene in gene_set:
				redundant_genes.append(shared_gene)
	
	return set(nonredundant_genes)

# Write output
for species in species_list:
	print("Writing output for %s..." % (species))
	for sample in pres_genes_dict['all'][species]:
		
		items = [species, sample, \
						len(pres_genes_dict['all'][species][sample]), \
						len(total_genes_dict['all'][species]), \
						len(sample_pres_genes_dict['all'][sample] - pres_genes_dict['all'][species][sample]), \
						len(sample_pres_gene_clusters_dict['all'][sample] - pres_gene_clusters_dict['all'][species][sample]), \
						len(pres_genes_dict['acc95'][species][sample]), \
						len(total_genes_dict['acc95'][species]), \
						len(sample_pres_genes_dict['acc95'][sample] - pres_genes_dict['acc95'][species][sample]), \
						len(sample_pres_gene_clusters_dict['acc95'][sample] - pres_gene_clusters_dict['acc95'][species][sample]), \
						len(pres_genes_dict['acc90'][species][sample]), \
						len(total_genes_dict['acc90'][species]), \
						len(sample_pres_genes_dict['acc90'][sample] - pres_genes_dict['acc90'][species][sample]), \
						len(sample_pres_gene_clusters_dict['acc90'][sample] - pres_gene_clusters_dict['acc90'][species][sample])]
		output_file.write(','.join(map(lambda x: str(x), items)) + '\n')

output_file.close()

# Dump pickles
pdir = '%s/pickles' % config.data_directory

# pickle.dump(gene_cluster_map, open('%s/gene_cluster_map.pkl' % pdir, 'w'))
# pickle.dump(pres_genes_dict, open('%s/pres_genes_dict.pkl' % pdir, 'wb'))
# pickle.dump(sample_pres_genes_dict, open('%s/sample_pres_genes_dict.pkl' % pdir, 'wb'))
# pickle.dump(total_genes_dict, open('%s/total_genes_dict.pkl' % pdir, 'wb'))
