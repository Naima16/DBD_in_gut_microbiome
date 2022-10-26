import config
import gzip
import pickle

pickle_dir = "%s/pickles" % config.data_directory

gene_gain_info_dict = pickle.load(open("%s/gene_gain_info_dict.pkl" % (pickle_dir), 'rb'))
gene_loss_info_dict = pickle.load(open("%s/gene_loss_info_dict.pkl" % (pickle_dir), 'rb'))

# Loss

genome_ids = set()
for sample_i, sample_j in gene_loss_info_dict:
	for species in gene_loss_info_dict[(sample_i, sample_j)]:
		losses = gene_loss_info_dict[(sample_i, sample_j)][species]
		for loss in losses:
			gene_name, D1, Dm1, D2, Dm2 = loss
			genome_id, gene = gene_name.split('.peg.')
			genome_ids.add(genome_id)

gene_id_desc_dict = {}
			
for genome_id in genome_ids:
	f = gzip.open('%s/features/%s.PATRIC.features.tab.gz' % (config.patric_directory, genome_id), 'r')
	header = f.readline() #header
	for line in f:
		items = line.strip().split("\t")
		if items[0] !='' and items[5] !='' and len(items)>14: # sometimes entries are blank
			gene_id =        items[5].split('|')[1] # id of gene
			gene_description = items[14] # what the gene does
			gene_id_desc_dict[gene_id] = gene_description # load into the dictionary

output_f = open("%s/poyet_gene_loss_details.tsv" % config.analysis_directory, 'w')
output_f.write('\t'.join(['species', 'sample1', 'sample2', 'gene_id', 'gene_desc', 'D1', 'Dm1', 'D2', 'Dm2']) + '\n')

for sample_i, sample_j in gene_loss_info_dict:
	for species in gene_loss_info_dict[(sample_i, sample_j)]:
		losses = gene_loss_info_dict[(sample_i, sample_j)][species]
		for loss in losses:
			gene_id, D1, Dm1, D2, Dm2 = loss
			gene_desc = gene_id_desc_dict[gene_id]
			output_f.write('\t'.join([str(val) for val in [species, sample_i, sample_j, gene_id, gene_desc, D1, Dm1, D2, Dm2]]) + '\n')

output_f.close()

# Gains

genome_ids = set()
for sample_i, sample_j in gene_gain_info_dict:
	for species in gene_gain_info_dict[(sample_i, sample_j)]:
		losses = gene_gain_info_dict[(sample_i, sample_j)][species]
		for loss in losses:
			gene_name, D1, Dm1, D2, Dm2 = loss
			genome_id, gene = gene_name.split('.peg.')
			genome_ids.add(genome_id)

gene_id_desc_dict = {}
			
for genome_id in genome_ids:
	f = gzip.open('%s/features/%s.PATRIC.features.tab.gz' % (config.patric_directory, genome_id), 'r')
	header = f.readline() #header
	for line in f:
		items = line.strip().split("\t")
		if items[0] !='' and items[5] !='' and len(items)>14: # sometimes entries are blank
			gene_id =        items[5].split('|')[1] # id of gene
			gene_description = items[14] # what the gene does
			gene_id_desc_dict[gene_id] = gene_description # load into the dictionary

output_f = open("%s/poyet_gene_gain_details.tsv" % config.analysis_directory, 'w')
output_f.write('\t'.join(['species', 'sample1', 'sample2', 'gene_id', 'gene_desc', 'D1', 'Dm1', 'D2', 'Dm2']) + '\n')

for sample_i, sample_j in gene_gain_info_dict:
	for species in gene_gain_info_dict[(sample_i, sample_j)]:
		losses = gene_gain_info_dict[(sample_i, sample_j)][species]
		for loss in losses:
			gene_id, D1, Dm1, D2, Dm2 = loss
			gene_desc = gene_id_desc_dict[gene_id]
			output_f.write('\t'.join([str(val) for val in [species, sample_i, sample_j, gene_id, gene_desc, D1, Dm1, D2, Dm2]]) + '\n')

output_f.close()
