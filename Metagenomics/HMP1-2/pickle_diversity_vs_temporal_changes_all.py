# Investigating 'amount of evolution' vs alpha diversity
# For HMP adults only
# Does not restrict to QP pairs
# This script pickles both SNP change and polymorphism info

from utils import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import sys, numpy as np, random, math
from utils import temporal_changes_utils
from collections import defaultdict
import bz2

rarefied_data_dir = "/u/home/d/daisyche/dbd/data_rarefied"

# ==========================================================
# Standard header to read in argument information
# ==========================================================
import argparse
parser = argparse.ArgumentParser()
# parser.add_argument('--species', type=str, help='Run the script for one specified species')
parser.add_argument('--sweep-type', type=str, help="Full or partial sweep", default="full")

args = parser.parse_args()
sweep_type = args.sweep_type

if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep-type. Choose from full, partial")

if sweep_type == 'full':
	lower_threshold, upper_threshold = 0.2, 0.8
elif sweep_type == 'partial':
	lower_threshold, upper_threshold = 0.35, 0.65
# ==========================================================

# Modification event threshold (config: 20)
modification_difference_threshold = config.modification_difference_threshold

# Replacement event threshold (config: 500)
replacement_difference_threshold = config.replacement_difference_threshold

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Samples
all_samples = su.get_sample_names('all')

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# Relative abundance file
relab_fpath = "%s/species/relative_abundance.txt.bz2" % config.data_directory
relab_file = bz2.BZ2File(relab_fpath, 'r')
data = [row.strip().split('\t') for row in relab_file]
samples = su.parse_merged_sample_names(data[0][1:])

# Generate alpha diversity, richness dictionary
alpha_div_dict = {}
richness_dict = {}
for i in range(len(samples)):
	acc, richness_acc = 0, 0
	for row in data[1:]:
		rel_ab = float(row[i+1]) # offset by 1
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
			richness_acc += 1
	alpha_div_dict[samples[i]] = (acc*-1)
	richness_dict[samples[i]] = richness_acc

# Relative abundance file (rarefied to 20 million)
relab_fpath = "%s/species/relative_abundance.txt" % rarefied_data_dir
relab_file = open(relab_fpath, 'r')
data = [row.strip().split('\t') for row in relab_file]
samples_rarefied = su.parse_merged_sample_names(data[0][1:])

# Generate alpha diversity, richness dictionary
alpha_div_dict_rarefied = {}
richness_dict_rarefied = {}
for i in range(len(samples_rarefied)):
	acc, richness_acc = 0, 0
	for row in data[1:]:
		rel_ab = float(row[i+1]) # offset by 1
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
			richness_acc += 1
	alpha_div_dict_rarefied[samples_rarefied[i]] = (acc*-1)
	richness_dict_rarefied[samples_rarefied[i]] = richness_acc

'''
# Get strain number info from Ricky's DBD folder
strain_num_dict = defaultdict(dict)
strain_num_cov10_dict = defaultdict(dict)

ricky_dbd_dir = '/u/project/ngarud/rwolff/dbd'
with open("%s/Poyet_strain_number.csv" % ricky_dbd_dir, 'r') as file:
	samples = file.readline().strip().split(',')[1:]
	for line in file:
		items = line.strip().split(',')
		species = items[0]
		for value, sample in zip(items[1:], samples):
			if value != '':
				strain_num_dict[sample][species] = int(value)

with open("%s/Poyet_strain_number_threshold_10.csv" % ricky_dbd_dir, 'r') as file:
	samples = file.readline().strip().split(',')[1:]
	for line in file:
		items = line.strip().split(',')
		species = items[0]
		for value, sample in zip(items[1:], samples):
			if value != '':
				strain_num_cov10_dict[sample][species] = int(value)
'''

# ====================================================================
# Dictionary: sample -> species -> within-sample polymorphism (average estimate)
sample_species_polymorphism_dict = defaultdict(dict)

# Dictionary: (sample1, sample2) pair -> species -> number of SNP changes
snp_change_dict = defaultdict(dict)

# Dictionary: (sample1, sample2) pair -> species -> gene gain/loss number
gene_gain_dict = defaultdict(dict)
gene_loss_dict = defaultdict(dict)

# Same as above but is entirely QP agnostic
all_gene_gain_dict = defaultdict(dict)
all_gene_loss_dict = defaultdict(dict)

# Dictionary: (sample1, sample2) pair -> species -> T/F QP sample pair
qp_pair_status_dict = defaultdict(dict)

# Dictionary: sample -> species -> T/F QP sample
qp_status_dict = defaultdict(dict)

# ====================================================================

for species_name in good_species_list:
	
	print("Working on " + species_name)
	
	# Store within-host polymorphism for ALL samples
	_, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D']))
	for sample in all_samples:		
		if sample not in sfs_map: # TODO
			continue
		within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
		try:
			within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05) # 95% CI
		except:
			continue
		sample_species_polymorphism_dict[sample][species_name] = (within_rate_lower + within_rate_upper)/2.0
	
	# Store temporal change info
	# Look at within-host QP sample pairs, choosing consecutive timepoint pairs
	# IN ADDITION, look at within-host sample pairs that are not QP
	qp_samples = sorted(su.calculate_qp_samples(all_samples, species_name)['qp'])
	
	for sample in all_samples:
		qp_status_dict[sample][species_name] = (sample in qp_samples)
	
	_, qp_same_subject_idxs, _ = su.calculate_ordered_subject_pairs(sample_order_map, qp_samples, within_host_type='consecutive')
	_, all_same_subject_idxs, _ = su.calculate_ordered_subject_pairs(sample_order_map, all_samples, within_host_type='consecutive')
	
	# Isolate sample pairs from all_same_subject_idxs where not both samples
	# are QP, so there isn't redundancy with qp_same_subject_idxs
	qp_all_idx_tuples = [] # QP sample pairs translated to all_samples idxs
	for qp_idx1, qp_idx2 in zip(qp_same_subject_idxs[0], qp_same_subject_idxs[1]):
		s1, s2 = qp_samples[qp_idx1], qp_samples[qp_idx2]
		qp_all_idx_tuples.append((all_samples.index(s1), all_samples.index(s2)))
	
	nonqp_same_subject_idxs = ([], []) # note that indices are w.r.t. all_samples
	for all_idx1, all_idx2 in zip(all_same_subject_idxs[0], all_same_subject_idxs[1]):
		if (all_idx1, all_idx2) not in qp_all_idx_tuples:
			nonqp_same_subject_idxs[0].append(all_idx1)
			nonqp_same_subject_idxs[1].append(all_idx2)
	
	# Load temporal change map
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name)
	
	# First look at SNP changes
	
	for desired_samples, same_subject_idxs, status in zip([qp_samples, all_samples, all_samples], [qp_same_subject_idxs, nonqp_same_subject_idxs, all_same_subject_idxs], ['qp', 'nonqp', 'all']):
		
		print("Looking at %i sample pairs..." % len(same_subject_idxs[0]))
		
		for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
			
			sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
			sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
			subject = sample_subject_map[sample_i]
			
			# Get SNP change info
			L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold, upper_threshold)
			
			nerr = L*perr
			num_mutations = len(mutations)
			num_reversions = len(reversions)
			num_snp_changes = num_mutations+num_reversions
			
			# Ignore iffy sample pairs
			print(L)
			if L < config.min_opportunities:
				print("Not enough opportunities")
				continue
			if (perr<-0.5) or (perr>0.5):
				print("Bad perr: %.02f" % perr)
				continue					
			if (nerr > max([0.5, 0.1*num_snp_changes])):
				print("Bad nerr: %.02f" % nerr)
				continue # Only take things with low-ish FPR
			
			# Get gene changes info
			gene_opps, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
			all_changes = gains + losses
			
			# Compute number of gene gains/losses
			num_gains = len(gains) if gains != None else 0
			num_losses = len(losses) if losses != None else 0
			num_gene_changes = num_gains + num_losses
			
			# Don't want to look at things with high error rates
			if ((gene_perr<-0.5) or (gene_perr>0.5)):
				continue
			
			# Store QP status
			if status == 'qp':
				qp_pair_status_dict[(sample_i, sample_j)][species_name] = True
			elif status == 'nonqp':
				qp_pair_status_dict[(sample_i, sample_j)][species_name] = False
			
			# Store gene change information
			gene_gain_dict[(sample_i, sample_j)][species_name] = num_gains
			gene_loss_dict[(sample_i, sample_j)][species_name] = num_losses
			
			if status == 'all':
				all_gene_gain_dict[(sample_i, sample_j)][species_name] = num_gains
				all_gene_loss_dict[(sample_i, sample_j)][species_name] = num_losses
			
			# Store number of SNP changes in modification event for host
			snp_change_dict[(sample_i, sample_j)][species_name] = num_snp_changes

# ====================================================================
# Pickle things
# ====================================================================

import pickle
pickle_dir = "%s/pickles" % config.data_directory

pickle.dump(all_gene_gain_dict, open("%s/all_gene_gain_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(all_gene_loss_dict, open("%s/all_gene_loss_dict.pkl" % (pickle_dir), 'wb'))
'''
pickle.dump(alpha_div_dict, open("%s/alpha_div_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(richness_dict, open("%s/richness_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(alpha_div_dict_rarefied, open("%s/alpha_div_dict_rarefied.pkl" % (pickle_dir), 'wb'))
pickle.dump(richness_dict_rarefied, open("%s/richness_dict_rarefied.pkl" % (pickle_dir), 'wb'))
pickle.dump(sample_species_polymorphism_dict, open("%s/sample_species_polymorphism_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(qp_pair_status_dict, open("%s/qp_pair_status_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(qp_status_dict, open("%s/qp_status_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(snp_change_dict, open("%s/snp_change_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(gene_gain_dict, open("%s/gene_gain_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(gene_loss_dict, open("%s/gene_loss_dict.pkl" % (pickle_dir), 'wb'))
# pickle.dump(strain_num_cov10_dict, open("%s/strain_num_cov10_dict.pkl" % (pickle_dir), 'wb'))
# pickle.dump(strain_num_dict, open("%s/strain_num_dict.pkl" % (pickle_dir), 'wb'))
'''

# ====================================================================
# Store table
# ====================================================================

analysis_dir = config.analysis_directory
with open("%s/HMP_temporal_changes_diversity_full.csv" % analysis_dir, 'w') as outfile:
	fields = ','.join('sample_tp1 sample_tp2 qp_status subject species num_SNP_changes num_gene_gains num_gene_losses alpha_div_tp1 alpha_div_tp2 richness_tp1 richness_tp2 polymorphism_tp1 polymorphism_tp2'.split())
	outfile.write(fields + '\n')
	for sample_i, sample_j in snp_change_dict:
		for species in snp_change_dict[(sample_i, sample_j)]:
			
			subject = sample_subject_map[sample_i]
			num_SNP_changes = snp_change_dict[(sample_i, sample_j)][species]
			num_gene_gains = gene_gain_dict[(sample_i, sample_j)][species]
			num_gene_losses = gene_loss_dict[(sample_i, sample_j)][species]
			qp_status = qp_pair_status_dict[(sample_i, sample_j)][species]
			adiv_tp1, adiv_tp2 = alpha_div_dict[sample_i], alpha_div_dict[sample_j]
			richness_tp1, richness_tp2 = richness_dict[sample_i], richness_dict[sample_j]
			polymorphism_tp1 = sample_species_polymorphism_dict[sample_i][species]
			polymorphism_tp2 = sample_species_polymorphism_dict[sample_j][species]
			
			fields_data = [sample_i, sample_j, qp_status, subject, species, num_SNP_changes, num_gene_gains, num_gene_losses, adiv_tp1, adiv_tp2, richness_tp1, richness_tp2, polymorphism_tp1, polymorphism_tp2]
			
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')
