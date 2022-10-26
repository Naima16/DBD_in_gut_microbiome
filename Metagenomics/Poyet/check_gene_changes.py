#!/usr/bin/env python
# coding: utf-8

# In[12]:


from utils import temporal_changes_utils, midas_db_utils, core_gene_utils, parse_patric, config
import pickle
from collections import defaultdict


# In[3]:


snp_change_counts = pickle.load(open('../data_poyet/pickles/snp_change_dict.pkl', 'rb'))
gene_gain_counts = pickle.load(open('../data_poyet/pickles/gene_gain_dict.pkl', 'rb'))

species_spair_gene_gain_dict = defaultdict(dict)

for tup in snp_change_counts:
	for species in gene_gain_counts[tup]:
		count = snp_change_counts[tup][species]
		if count <= 20 and gene_gain_counts[tup][species] > 0:
			species_spair_gene_gain_dict[species][tup] = gene_gain_counts[tup][species]


# In[ ]:


species_spair_gene_gain_full_dict = defaultdict(dict)

for species in species_spair_gene_gain_dict:
	print(species)
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species)
	for sample_i, sample_j in species_spair_gene_gain_dict[species]:
		gene_opps, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
		if len(gains) != species_spair_gene_gain_dict[species][(sample_i, sample_j)]:
			print("Huh?")
		species_spair_gene_gain_full_dict[species][(sample_i, sample_j)] = gains


# In[5]:


print(species_spair_gene_gain_full_dict)


# In[8]:


species_spair_gene_gain_names_dict = defaultdict(dict)

for species in species_spair_gene_gain_full_dict:
    print(species)
    # get all genome ids for this species' pan genome:
    genome_ids = midas_db_utils.get_ref_genome_ids(species)
    
    # Load the non-shared genes (whitelisted genes):
    non_shared_genes = core_gene_utils.parse_non_shared_pangenome_genes(species)

    # load the gene descriptions for all genomes coresponding to this speceis:
    gene_descriptions=parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
    
    for spair in species_spair_gene_gain_full_dict[species]:
        gene_names = []
        for gene_id, D1, Dm1, D2, Dm2 in species_spair_gene_gain_full_dict[species][spair]:
            try:
                gene_names.append(gene_descriptions[gene_id])
            except:
                print(gene_id)
        species_spair_gene_gain_names_dict[species][spair] = gene_names


# In[19]:


f = open('%s/gene_gain_annotations.tsv' % config.analysis_directory, 'w')
f.write('\t'.join(['sample1', 'sample2', 'species', 'gain_count', 'gain_genes']) + '\n')
for species in species_spair_gene_gain_names_dict:
    for spair in species_spair_gene_gain_names_dict[species]:
        gain_count = species_spair_gene_gain_dict[species][spair]
        gain_genes = species_spair_gene_gain_names_dict[species][spair]
        f.write('\t'.join([spair[0], spair[1], species, str(gain_count), ','.join(gain_genes)]) + '\n')
f.close()


# In[9]:


snp_change_counts = pickle.load(open('../data_poyet/pickles/snp_change_dict.pkl', 'rb'))
gene_loss_counts = pickle.load(open('../data_poyet/pickles/gene_loss_dict.pkl', 'rb'))

species_spair_gene_loss_dict = defaultdict(dict)

for tup in snp_change_counts:
	for species in gene_loss_counts[tup]:
		count = snp_change_counts[tup][species]
		if count <= 20 and gene_loss_counts[tup][species] > 0:
			species_spair_gene_loss_dict[species][tup] = gene_loss_counts[tup][species]


# In[ ]:


species_spair_gene_loss_full_dict = defaultdict(dict)

for species in species_spair_gene_loss_dict:
	print(species)
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species)
	for sample_i, sample_j in species_spair_gene_loss_dict[species]:
		gene_opps, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
		if len(losses) != species_spair_gene_loss_dict[species][(sample_i, sample_j)]:
			print("Huh?")
		species_spair_gene_loss_full_dict[species][(sample_i, sample_j)] = losses


# In[15]:


print(species_spair_gene_loss_full_dict)


# In[16]:


species_spair_gene_loss_names_dict = defaultdict(dict)

for species in species_spair_gene_loss_full_dict:
    print(species)
    # get all genome ids for this species' pan genome:
    genome_ids = midas_db_utils.get_ref_genome_ids(species)
    
    # Load the non-shared genes (whitelisted genes):
    non_shared_genes = core_gene_utils.parse_non_shared_pangenome_genes(species)

    # load the gene descriptions for all genomes coresponding to this speceis:
    gene_descriptions=parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
    
    for spair in species_spair_gene_loss_full_dict[species]:
        gene_names = []
        for gene_id, D1, Dm1, D2, Dm2 in species_spair_gene_loss_full_dict[species][spair]:
            try:
                gene_names.append(gene_descriptions[gene_id])
            except:
                print(gene_id)
        species_spair_gene_loss_names_dict[species][spair] = gene_names


# In[20]:


f = open('%s/gene_loss_annotations.tsv' % config.analysis_directory, 'w')
f.write('\t'.join(['sample1', 'sample2', 'species', 'loss_count', 'loss_genes']) + '\n')
for species in species_spair_gene_loss_names_dict:
    for spair in species_spair_gene_loss_names_dict[species]:
        loss_count = species_spair_gene_loss_dict[species][spair]
        loss_genes = species_spair_gene_loss_names_dict[species][spair]
        f.write('\t'.join([spair[0], spair[1], species, str(loss_count), ','.join(loss_genes)]) + '\n')
f.close()

