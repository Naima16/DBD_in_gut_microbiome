from utils import config, sample_utils as su, parse_midas_data, temporal_changes_utils, substitution_rates_utils
from numpy.random import choice
import sys, pickle
from collections import defaultdict

# Parameters ==================================================================
lower_threshold, upper_threshold = 0.2, 0.8

min_sample_size = 3
variant_types = ['1D','4D']
within_host_type = 'nonconsecutive' # nonconsecutive timepoints
min_snp_change_sample_size = 5

modification_difference_threshold = config.modification_difference_threshold
replacement_difference_threshold = config.replacement_difference_threshold
# =============================================================================

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_subject_map = su.parse_sample_subject_map()
sample_order_map = su.parse_sample_order_map()
sys.stderr.write("Done!\n")

# HMP samples
hmp_samples = sample_subject_map.keys()

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# Pooled SNP/gene change distributions
pooled_snp_change_distribution = defaultdict(list)
pooled_gene_change_distribution = defaultdict(list)
pooled_between_snp_change_distribution = defaultdict(list)
pooled_between_gene_change_distribution = defaultdict(list)

# Pickle file names
pooled_snp_pickle_fn = '%s/pickles/pooled_snp_change.pkl' % (config.data_directory)
pooled_gene_pickle_fn = '%s/pickles/pooled_gene_change.pkl' % (config.data_directory)
pooled_between_snp_pickle_fn = '%s/pickles/pooled_between_snp_change.pkl' % (config.data_directory)
pooled_between_gene_pickle_fn = '%s/pickles/pooled_between_gene_change.pkl' % (config.data_directory)

for species_name in good_species_list:
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	qp_sample_set = su.load_qp_samples(hmp_samples, species_name)['qp']
	qp_sample_list = list(qp_sample_set)
	
	# Must have minimum of 3 QP samples to proceed
	if len(qp_sample_list) < min_sample_size:
		continue
	
	# Load substitution rates
	substitution_rate_map = substitution_rates_utils.load_substitution_rate_map(species_name)
	
	if substitution_rate_map == {}: # Not enough haploid samples
		sys.stderr.write("Not enough haploid samples!\n")
		continue
	
	sys.stderr.write("Calculating SNV matrix...\n")
	dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=qp_sample_list)
	snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
	snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
	snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
	sys.stderr.write("Done!\n")
	
	sys.stderr.write("Loading gene matrix...\n")
	gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=qp_sample_list)
	gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
	gene_opportunity_matrix = gene_loss_opportunity_matrix
	gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
	sys.stderr.write("Done!\n")
	
	# Load temporal changes
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name)
	sys.stderr.write("Done!\n")
	
	desired_samples = qp_sample_list
	same_subject_idxs = su.calculate_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type)
	
	# Loop over different pairs of within-host samples
	for idx_i, idx_j in zip(same_subject_idxs[0], same_subject_idxs[1]):			 
		
		sample_i, sample_j = desired_samples[idx_i], desired_samples[idx_j]
		
		tp_i = 'A' + str(sample_order_map[sample_i][1])
		tp_j = 'A' + str(sample_order_map[sample_j][1])
		tp_pair = frozenset((tp_i, tp_j))
		
		good_idxs = su.calculate_samples_in_different_subjects(sample_subject_map, qp_sample_list, sample_i)
		good_idxs *= ((snp_opportunity_matrix[idx_i,:]>0.5) * (gene_opportunity_matrix[idx_j,:]>0.5))
		
		if good_idxs.sum() < 1:
			# print("BAM (%i, %i)" % (idx_i, idx_j))
			continue
		
		L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
		
		if L<config.min_opportunities:
			# print("KERFLUFFLE (%i, %i)" % (idx_i, idx_j))
			continue
		
		nerr = L*perr
		
		num_mutations = len(mutations)
		num_reversions = len(reversions)
		num_snp_changes = num_mutations+num_reversions
		
		gene_L, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) # min_normal_copynum = 0.6, max_normal_copynum = 1.2)
		
		gene_nerr = gene_L*gene_perr
		num_gains = len(gains)
		num_losses = len(losses)
		num_gene_changes = num_gains+num_losses
		
		if (perr<-0.5) or (gene_perr < -0.5):
			# print("WHAT THE HECK (%i, %i)" % (idx_i, idx_j))
			continue
		
		if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
			# print("AIYA (%i, %i)" % (idx_i, idx_j))
			continue # Only take things with low-ish FPR
		
		# Pooled distributions
		pooled_snp_change_distribution[tp_pair].append(num_snp_changes)
		pooled_gene_change_distribution[tp_pair].append(num_gene_changes)
		pooled_between_snp_change_distribution[tp_pair].append(choice(snp_difference_matrix[idx_i, good_idxs])) 
		pooled_between_gene_change_distribution[tp_pair].append(choice(gene_difference_matrix[idx_i, good_idxs]))

# Dump the pickles!
pickle.dump(pooled_snp_change_distribution, open(pooled_snp_pickle_fn, 'wb'))
pickle.dump(pooled_gene_change_distribution, open(pooled_gene_pickle_fn, 'wb'))
pickle.dump(pooled_between_snp_change_distribution, open(pooled_between_snp_pickle_fn, 'wb'))
pickle.dump(pooled_between_gene_change_distribution, open(pooled_between_gene_pickle_fn, 'wb'))