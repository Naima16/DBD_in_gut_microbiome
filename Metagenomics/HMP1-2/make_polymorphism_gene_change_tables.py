# Investigating 'amount of evolution' vs alpha diversity
# For HMP adults only
# Does not restrict to QP pairs
# Construct final tables from pickles from pickle_diversity_vs_temporal_changes_all.py

from utils import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import sys, numpy as np, random, math
from utils import temporal_changes_utils
from collections import defaultdict
import bz2

# ====================================================================
# Load pickles etc.
# ====================================================================

import pickle
pickle_dir = "%s/pickles" % config.data_directory

sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()

alpha_div_dict = pickle.load(open("%s/alpha_div_dict.pkl" % (pickle_dir), 'rb'))
richness_dict = pickle.load(open("%s/richness_dict.pkl" % (pickle_dir), 'rb'))
alpha_div_dict_rarefied = pickle.load(open("%s/alpha_div_dict_rarefied.pkl" % (pickle_dir), 'rb'))
richness_dict_rarefied = pickle.load(open("%s/richness_dict_rarefied.pkl" % (pickle_dir), 'rb'))
sample_species_polymorphism_dict = pickle.load(open("%s/sample_species_polymorphism_dict.pkl" % (pickle_dir), 'rb'))
qp_pair_status_dict = pickle.load(open("%s/qp_pair_status_dict.pkl" % (pickle_dir), 'rb'))
qp_status_dict = pickle.load(open("%s/qp_status_dict.pkl" % (pickle_dir), 'rb'))
snp_change_dict = pickle.load(open("%s/snp_change_dict.pkl" % (pickle_dir), 'rb'))
gene_gain_dict = pickle.load(open("%s/all_gene_gain_dict.pkl" % (pickle_dir), 'rb'))
gene_loss_dict = pickle.load(open("%s/all_gene_loss_dict.pkl" % (pickle_dir), 'rb'))

# ====================================================================
# Polymorphism per sample-species pair table (HMP)
# Also includes alpha diversity and richness
# ====================================================================

analysis_dir = config.analysis_directory
with open("%s/HMP_polymorphism.csv" % analysis_dir, 'w') as outfile:
	
	fields = ','.join('sample_id subject tp_order species alpha_div alpha_div_rare richness richness_rare polymorphism'.split())
	outfile.write(fields + '\n')
	
	for sample in sample_species_polymorphism_dict:
		for species in sample_species_polymorphism_dict[sample]:
			
			subject, order = sample_order_map[sample]
			qp_status = qp_status_dict[sample][species]
			alpha_div = alpha_div_dict[sample]
			alpha_div_rarefied = alpha_div_dict_rarefied[sample]
			richness = richness_dict[sample]
			richness_rarefied = richness_dict_rarefied[sample]
			polymorphism = sample_species_polymorphism_dict[sample][species]
			
			fields_data = [sample, subject, order, species, alpha_div, alpha_div_rarefied, richness, richness_rarefied, polymorphism]
			
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')

# ====================================================================
# Delta polymorphism for consecutive timepoint pairs (HMP)
# ====================================================================

desired_samples = sample_species_polymorphism_dict.keys()
same_subject_idxs = su.calculate_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type='consecutive')

analysis_dir = config.analysis_directory
with open("%s/HMP_delta_polymorphism.csv" % analysis_dir, 'w') as outfile:
	
	fields = ','.join('sample_tp1 sample_tp2 subject tp1 tp2 species delta_polymorphism polymorphism_tp1 polymorphism_tp2 alpha_div_tp1 alpha_div_tp2 alpha_div_rare_tp1 alpha_div_rare_tp2 richness_tp1 richness_tp2 richness_rare_tp1 richness_rare_tp2 '.split())
	outfile.write(fields + '\n')
	
	for i, j in zip(same_subject_idxs[0], same_subject_idxs[1]):
		sample_i = desired_samples[i]
		sample_j = desired_samples[j]
		
		shared_species = set(sample_species_polymorphism_dict[sample_i].keys()).intersection(set(sample_species_polymorphism_dict[sample_j].keys()))
		
		for species in shared_species:
			
			subject = sample_subject_map[sample_i]
			subject1, tp1 = sample_order_map[sample_i]
			subject2, tp2 = sample_order_map[sample_j]
			if subject1 != subject2 or subject != subject1:
				print("Weird")
			adiv_tp1, adiv_tp2 = alpha_div_dict[sample_i], alpha_div_dict[sample_j]
			adiv_rare_tp1, adiv_rare_tp2 = alpha_div_dict_rarefied[sample_i], alpha_div_dict_rarefied[sample_j]
			richness_tp1, richness_tp2 = richness_dict[sample_i], richness_dict[sample_j]
			richness_rare_tp1, richness_rare_tp2 = richness_dict_rarefied[sample_i], richness_dict_rarefied[sample_j]
			polymorphism_tp1 = sample_species_polymorphism_dict[sample_i][species]
			polymorphism_tp2 = sample_species_polymorphism_dict[sample_j][species]
			delta_polymorphism = polymorphism_tp2 - polymorphism_tp1
			
			fields_data = [sample_i, sample_j, subject, tp1, tp2, species, delta_polymorphism, polymorphism_tp1, polymorphism_tp2, adiv_tp1, adiv_tp2, adiv_rare_tp1, adiv_rare_tp2, richness_tp1, richness_tp2, richness_rare_tp1, richness_rare_tp2]
			
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')

# ====================================================================
# Gene changes per sample pair (HMP)
# ====================================================================

with open("%s/HMP_gene_changes_full.csv" % analysis_dir, 'w') as outfile:
	fields = ','.join('sample_tp1 sample_tp2 subject tp1 tp2 species num_gene_gains num_gene_losses alpha_div_tp1 alpha_div_tp2 alpha_div_rare_tp1 alpha_div_rare_tp2 richness_tp1 richness_tp2 richness_rare_tp1 richness_rare_tp2 polymorphism_tp1 polymorphism_tp2'.split())
	outfile.write(fields + '\n')
	for sample_i, sample_j in gene_gain_dict:
		for species in gene_gain_dict[(sample_i, sample_j)]:
			
			subject = sample_subject_map[sample_i]
			subject1, tp1 = sample_order_map[sample_i]
			subject2, tp2 = sample_order_map[sample_j]
			if subject1 != subject2 or subject != subject1:
				print("Weird")
			num_SNP_changes = snp_change_dict[(sample_i, sample_j)][species]
			num_gene_gains = gene_gain_dict[(sample_i, sample_j)][species]
			num_gene_losses = gene_loss_dict[(sample_i, sample_j)][species]
			# qp_status = qp_pair_status_dict[(sample_i, sample_j)][species]
			adiv_tp1, adiv_tp2 = alpha_div_dict[sample_i], alpha_div_dict[sample_j]
			adiv_rare_tp1, adiv_rare_tp2 = alpha_div_dict_rarefied[sample_i], alpha_div_dict_rarefied[sample_j]
			richness_tp1, richness_tp2 = richness_dict[sample_i], richness_dict[sample_j]
			richness_rare_tp1, richness_rare_tp2 = richness_dict_rarefied[sample_i], richness_dict_rarefied[sample_j]
			polymorphism_tp1 = sample_species_polymorphism_dict[sample_i][species]
			polymorphism_tp2 = sample_species_polymorphism_dict[sample_j][species]
			
			fields_data = [sample_i, sample_j, subject, tp1, tp2, species, num_gene_gains, num_gene_losses, adiv_tp1, adiv_tp2, adiv_rare_tp1, adiv_rare_tp2, richness_tp1, richness_tp2, richness_rare_tp1, richness_rare_tp2, polymorphism_tp1, polymorphism_tp2]
			
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')
