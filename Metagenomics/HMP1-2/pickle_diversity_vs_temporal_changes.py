# Investigating 'amount of evolution' vs alpha diversity
# For HMP adults only
# This script pickles both SNP change and polymorphism info

from utils import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import pylab, sys, numpy as np, random, math
from utils import temporal_changes_utils
from collections import defaultdict
import bz2

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
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Samples
hmp_samples = su.get_sample_names('all') # Note: c's removed

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# Relative abundance file
relab_fpath = "%s/species/relative_abundance.txt.bz2" % config.data_directory
relab_file = bz2.BZ2File(relab_fpath, 'r')
data = [row.strip().split('\t') for row in relab_file]
samples = su.parse_merged_sample_names(data[0][1:])

# Generate alpha diversity dictionary
alpha_div_dict = {}
for i in range(len(samples)):
	acc = 0
	for row in data[1:]:
		rel_ab = float(row[i+1]) # offset by 1
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
	alpha_div_dict[samples[i]] = (acc*-1)

# ====================================================================
# Dictionary: host -> number of SNP changes aggregated across species
# Only include if modification event happened in host (longest)
host_snp_change_dict = defaultdict(int)

# Dictionary: host-species pair -> number of SNP changes
host_species_snp_change_dict = defaultdict(int)

# Dictionary: hist-species pair -> sample pair (longest)
host_species_sample_pair_dict = {}

# Dictionary: host -> alpha diversity at first timepoint
# Only include if modification event happened in host (longest)
host_alpha_diversity_dict = {}

# Dictionary: none / mod / replace -> host -> number of events [i.e. species with the relevant event type, for longest timepoint pair]
host_change_type_dict = {type: defaultdict(int) for type in ['none', 'mod', 'replace']}

# Dictionary: sample(-species pair) -> within-sample polymorphism (average estimate)
sample_species_polymorphism_dict = defaultdict(dict)

# Dictionary: host -> gene gain/loss numbers aggregated across species
host_gene_gain_dict = defaultdict(int)
host_gene_loss_dict = defaultdict(int)

# Dictionary: host-species pair -> gene gain/loss numbers
host_species_gene_gain_dict = defaultdict(int)
host_species_gene_loss_dict = defaultdict(int)
# ====================================================================

for species_name in good_species_list:
	
	print("Working on " + species_name)
	
	# Store within-host polymorphism for ALL samples
	_, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D']))
	for sample in hmp_samples:		
		if sample not in sfs_map: # TODO
			continue
		within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
		try:
			within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05) # 95% CI
		except:
			continue
		sample_species_polymorphism_dict[sample][species_name] = (within_rate_lower + within_rate_upper)/2.0
	
	# Store temporal change info
	# Only looking at within-host QP sample pairs, choosing longest
	# timepoint pair when there are multiple
	desired_samples = sorted(su.calculate_qp_samples(hmp_samples, species_name)['qp'])
	_, same_subject_idxs, _ = su.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type='longest')
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name)
	
	# First look at SNP changes
	for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
		
		sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
		sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
		
		subject = sample_subject_map[sample_i]
		
		L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold, upper_threshold)
		
		nerr = L*perr
		num_mutations = len(mutations)
		num_reversions = len(reversions)
		num_snp_changes = num_mutations+num_reversions
		
		# Ignore iffy sample pairs
		if L < config.min_opportunities:
			continue
		if (perr<-0.5) or (perr>0.5):
			continue					
		if (nerr > max([0.5, 0.1*num_snp_changes])):
			continue # Only take things with low-ish FPR
		if num_snp_changes == 0:
			host_change_type_dict['none'][subject] += 1
		elif num_snp_changes >= replacement_difference_threshold:
			host_change_type_dict['replace'][subject] += 1
		else:
			# Count anything not replacement and nonzero as modification here
			host_change_type_dict['mod'][subject] += 1
		
		# Store alpha (Shannon) diversity for earlier sample
		host_alpha_diversity_dict[subject] = (alpha_div_dict[sample_i], alpha_div_dict[sample_j])
		
		host_species_sample_pair_dict[(subject, species_name)] = (sample_i, sample_j)
		
		# Look at gene changes for non-replacement events
		if num_snp_changes < replacement_difference_threshold:
			
			gene_opps, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
			all_changes = gains + losses
			
			# Compute number of gene gains/losses
			num_gains = len(gains) if gains != None else 0
			num_losses = len(losses) if losses != None else 0
			num_gene_changes = num_gains + num_losses
			
			# Don't want to look at things with high error rates
			if not ((gene_perr<-0.5) or (gene_perr>0.5)):
				# Store gene change information
				host_gene_gain_dict[subject] += num_gains
				host_gene_loss_dict[subject] += num_losses
				host_species_gene_gain_dict[(subject, species_name)] = num_gains
				host_species_gene_loss_dict[(subject, species_name)] = num_losses
		
		# Store number of SNP changes in modification event for host
		if num_snp_changes < modification_difference_threshold:
			host_snp_change_dict[subject] += num_snp_changes
			host_species_snp_change_dict[(subject, species_name)] = num_snp_changes
		
		# Store alpha (Shannon) diversity for earlier sample
		host_alpha_diversity_dict[subject] = (alpha_div_dict[sample_i], alpha_div_dict[sample_j])

# Pickle things

import pickle
pickle_dir = "%s/pickles" % config.data_directory

pickle.dump(alpha_div_dict, open("%s/alpha_div_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_snp_change_dict, open("%s/host_snp_change_dict_%s.pkl" % (pickle_dir, sweep_type), 'wb'))
pickle.dump(host_species_snp_change_dict, open("%s/host_species_snp_change_dict_%s.pkl" % (pickle_dir, sweep_type), 'wb'))
pickle.dump(host_alpha_diversity_dict, open("%s/host_alpha_diversity_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_change_type_dict, open("%s/host_change_type_dict_%s.pkl" % (pickle_dir, sweep_type), 'wb'))
pickle.dump(sample_species_polymorphism_dict, open("%s/sample_species_polymorphism_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_gene_gain_dict, open("%s/host_gene_gain_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_gene_loss_dict, open("%s/host_gene_loss_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_species_gene_gain_dict, open("%s/host_species_gene_gain_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_species_gene_loss_dict, open("%s/host_species_gene_loss_dict.pkl" % (pickle_dir), 'wb'))
pickle.dump(host_species_sample_pair_dict, open("%s/host_species_sample_pair_dict.pkl" % (pickle_dir), 'wb'))

# Store polymorphism and SNP change tables

analysis_dir = config.analysis_directory

with open("%s/HMP1-2_polymorphism_rates_alpha_divs.csv" % analysis_dir, 'w') as outfile:
	fields = ','.join('sample_id subject_id species_name polymorphism_rate Shannon_alpha_diversity'.split())
	outfile.write(fields + '\n')
	for sample in sample_species_polymorphism_dict:
		subject = sample_subject_map[sample]
		alpha_div = alpha_div_dict[sample]
		for species in sample_species_polymorphism_dict[sample]:
			prate = sample_species_polymorphism_dict[sample][species]
			fields_data = [sample, subject, species, prate, alpha_div]
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')

with open("%s/HMP1-2_SNP_changes_alpha_divs.csv" % analysis_dir, 'w') as outfile:
	fields = ','.join('subject_id species_name num_SNP_changes sample_id_tp1 sample_id_tp2 alpha_div_tp1 alpha_div_tp2'.split())
	outfile.write(fields + '\n')
	for subject, species in host_species_snp_change_dict:
		num_SNP_changes = host_species_snp_change_dict[(subject, species)]
		sample_i, sample_j = host_species_sample_pair_dict[(subject, species)]
		adiv_tp1, adiv_tp2 = host_alpha_diversity_dict[subject]
		fields_data = [subject, species, num_SNP_changes, sample_i, sample_j, adiv_tp1, adiv_tp2]
		outfile.write(','.join([str(f) for f in fields_data]) + '\n')

with open('%s/HMP1-2_gene_changes_alpha_divs.csv' % analysis_dir, 'w') as outfile:
	outfile.write('subject_id,species_id,num_gene_gains,num_gene_losses,sample_id_tp1,sample_id_tp2,alpha_div_tp1,alpha_div_tp2\n')
	for host, species in host_species_gene_gain_dict:
		sample_i, sample_j = host_species_sample_pair_dict[(host, species)]
		gene_gains = host_species_gene_gain_dict[(host, species)]
		gene_losses = host_species_gene_loss_dict[(host, species)]
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[host]		
		outfile.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (host, species, gene_gains, gene_losses, sample_i, sample_j, alpha_div1, alpha_div2))
