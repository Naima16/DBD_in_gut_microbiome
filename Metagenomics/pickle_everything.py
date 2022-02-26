from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils, core_gene_utils
import numpy as np
from numpy.random import choice, random, randint
from collections import defaultdict
import pickle
import sys

# ======================================================
# Examines all consecutive timepoint pairs within hosts
# across all cohorts, and pickles SNP/gene change info
# ======================================================

# Parameters
sweep_type = 'full' # assume full for now
thresholds = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
lower_threshold, upper_threshold = thresholds[sweep_type]

min_sample_size = 3
min_haploid_sample_size = 10
variant_types = ['1D','4D']
within_host_type = 'consecutive' # consecutive timepoints
min_snp_change_sample_size = 5

# For partitioning SNVs according to prevalence
derived_freq_bins = np.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = np.arange(0,len(derived_freq_bins)-1)
derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']

# For partitioning genes into different prevalence classes
gene_freq_bins = np.array([-1,0.1,0.5,0.9,2])
gene_freq_xticks			= [-4, -3,	-2,		-1,		0,	 1,		 2,		3, 4]
gene_freq_xticklabels = ['0','0.1','0.5', '0.9','1','0.9','0.5', '0.1','0']
gene_gain_virtual_freqs = np.array([3.5,2.5,1.5,0.5])
gene_loss_virtual_freqs = np.array([-3.5,-2.5,-1.5,-0.5])

# Sample-subject-order maps
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Cohorts
cohorts = ['hmp']

# Prevalence cohorts
prev_cohorts = ['all']

# Samples for each cohort
samples = {'hmp': su.get_sample_names()}

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

def get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map):
 
		gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
								
		f1 = A1*1.0/D2
		f2 = A2*1.0/D2
								
		is_reversion = (f1>f2)
								
		location_tuple = (contig, position)
								
		is_private_snv = (location_tuple in private_snv_map)
								
										
		# Now calculate frequency-stratified version
								
		if location_tuple in snv_freq_map:
				f = snv_freq_map[location_tuple]
		else:
				sys.stderr.write("SNP not in map. Shouldn't happen!\n")
				f = -0.5
								
		# Let's impose that private snvs have zero freq (specifically, lim 0^-)				 
		if is_private_snv:
				f = -0.5
								
		# Change f so that it represents
		# frequency of allele at second timepoint
		if is_reversion:
				f = 1-f
				
		return f

# ===================================================================
# Species SNP/gene change distributions

# species -> (sample1, sample2) -> [list of SNP change tuples]
# OR # SNP changes if that number exceeds 20 (replacement)
# where each SNP change tuple consists of
# (gene_name, contig, position, variant_type, A1, D1, A2, D2)
snp_changes = {species: {} for species in good_species_list}

# species -> (sample1, sample2) -> (gene gain tuples, gene loss tuples)
# OR (# gene gains, # gene losses) if replacement event
# where each gene change tuple consists of...
gene_changes = {species: {} for species in good_species_list}

# species -> (sample1, sample2) -> [list of (variant_type, prevalence, opps) tuples]
# only for modification events
snp_change_freqs = {species: defaultdict(list) for species in good_species_list}
snp_change_null_freqs = {species: defaultdict(list) for species in good_species_list}

# species -> (sample1, sample2) -> [list of prevalences]
# only for modification events
gene_gain_freqs = {species: defaultdict(list) for species in good_species_list}
gene_loss_freqs = {species: defaultdict(list) for species in good_species_list}
gene_loss_null_freqs = {species: defaultdict(list) for species in good_species_list}

# ===================================================================

for species_name in good_species_list[::-1]:
	
	# just for testing
	if species_name in ['Escherichia_coli_58110', 'Bacteroides_vulgatus_57955', 'Enterococcus_faecalis_56297']:
		continue
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	# Grab QP samples for this species
	qp_sample_lists = {}
	for cohort in cohorts:
		qp_sample_lists[cohort] = sorted(su.load_qp_samples(samples[cohort], species_name)['qp'])
	
	combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
	combined_sample_idx_map = {combined_qp_samples[i] : i for i in range(len(combined_qp_samples))}
	
	# Using all QP samples to threshold on sample size
	if len(combined_qp_samples) < min_haploid_sample_size: # 10
		sys.stderr.write("Not enough haploid samples!\n")
		continue
	
	# Load substitution rates for all QP samples
	sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
	substitution_rate_map = substitution_rates_utils.load_substitution_rate_map(species_name)
	
	if substitution_rate_map == {}: # Not enough haploid samples
		sys.stderr.write("Not enough haploid samples!\n")
		continue
	
	sys.stderr.write("Calculating SNV matrix...\n")
	dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_qp_samples)
	
	snp_difference_matrix = snp_mut_difference_matrix + snp_rev_difference_matrix
	snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
	snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
	sys.stderr.write("Done!\n")
	
	sys.stderr.write("Loading gene matrix...\n")
	gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_qp_samples)
	gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
	gene_opportunity_matrix = gene_loss_opportunity_matrix
	gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
	sys.stderr.write("Done!\n")
	
	sys.stderr.write("Loading 1D & 4D opportunity matrices...\n")	
	difference_matrices, opportunity_matrices = {}, {}	
	for var_type in variant_types:		
		matrix_samples, difference_matrix, opportunity_matrix = substitution_rates_utils.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=combined_qp_samples)		
		difference_matrices[var_type] = difference_matrix
		opportunity_matrices[var_type] = opportunity_matrix
	
	sys.stderr.write("Done!\n")
	
	# Load temporal change map
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name)
	
	# Load prevalences
	snv_freq_map = {prev_cohort: snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=True) for prev_cohort in prev_cohorts}
	
	# Load gene frequencies
	gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
	gene_freq_values = np.array(gene_freq_map.values())
	gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()
	
	# Loop over different cohorts
	for cohort in cohorts:		
		desired_samples = qp_sample_lists[cohort]
		
		same_subject_idxs = su.calculate_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type)
		
		# Get number of QP samples involved in same host pairs
		# for specific cohort, to replicate filter in old scripts
		# Loop over different pairs of within-host samples
		
		num_qp_within_host_samples = set()
		
		for sample_pair_idx in range(len(same_subject_idxs[0])):			 
			num_qp_within_host_samples.add(desired_samples[same_subject_idxs[0][sample_pair_idx]])
			num_qp_within_host_samples.add(desired_samples[same_subject_idxs[1][sample_pair_idx]])
		
		if len(num_qp_within_host_samples) < min_sample_size:
			continue
		
		# Loop over different pairs of within-host samples
		for sample_pair_idx in range(len(same_subject_idxs[0])):			 
				
				sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
				sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
				subject = (sample_subject_map[sample_i], sample_subject_map[sample_j])
				
				i = combined_sample_idx_map[sample_i]
				j = combined_sample_idx_map[sample_j]
				
				# Checks if among those samples from different hosts,
				# at least one of them has nonzero SNP and gene opportunities
				good_idxs = su.calculate_samples_in_different_subjects(sample_subject_map, combined_qp_samples, sample_i)
				good_idxs *= ( (snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5) )
				
				# FIRST FILTER
				if good_idxs.sum() < 1:
					sys.stderr.write("Not enough other-host samples!\n")
					continue
				
				matrix_idx_i = matrix_samples.index(sample_i)
				matrix_idx_j = matrix_samples.index(sample_j)
				
				# Numbers of site differences and opportunities between the timepoints
				nonsyn_diffs = difference_matrices['1D'][matrix_idx_i][matrix_idx_j]
				nonsyn_opps = opportunity_matrices['1D'][matrix_idx_i][matrix_idx_j]		
				syn_diffs = difference_matrices['4D'][matrix_idx_i][matrix_idx_j]
				syn_opps = opportunity_matrices['4D'][matrix_idx_i][matrix_idx_j]
				
				# SNP temporal changes
				L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
				
				# SECOND FILTER
				if L<config.min_opportunities:
					sys.stderr.write("Not enough SNP opportunities (should be >=100,000)!\n")
					continue
				
				nerr = L*perr
				
				num_mutations = len(mutations)
				num_reversions = len(reversions)
				num_snp_changes = num_mutations + num_reversions
				
				# Gene temporal changes
				gene_L, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
				
				gene_nerr = gene_L*gene_perr
				num_gains = len(gains)
				num_losses = len(losses)
				num_gene_changes = num_gains+num_losses
				
				# THIRD FILTER
				if (perr<-0.5) or (gene_perr < -0.5):
					sys.stderr.write("Perr too high!\n")
					continue
				
				# FOURTH FILTER
				if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
					sys.stderr.write("Nerr too high!\n")
					continue # Only take things with low-ish FPR
				
				# Store information
				if num_snp_changes <= 20: # Modification event
					
					snp_changes[species_name][(sample_i, sample_j)] = (mutations + reversions) #!
					gene_changes[species_name][(sample_i, sample_j)] = (gains, losses) #!
					
					# Delve into prevalence of modified SNPs and genes
					null_freq_dict = defaultdict(list)
					
					for snp_change in (mutations + reversions):
						
						variant_type = snp_change[3]
						freq_dict = {}
						for prev_cohort in prev_cohorts:
							f = get_sweep_prevalence(snp_change, snv_freq_map[prev_cohort], private_snv_map)
							freq_dict[prev_cohort] = f
						
						# Calculate null version based on # of opportunities
						opp_dict = {}
						total_opportunities = sum([opportunity_matrices[vt][i,j] for vt in variant_types])
						
						for vt in variant_types:
							opp_dict[vt] = opportunity_matrices[vt][i,j]/total_opportunities
						
						syn_opps = opportunity_matrices['4D'][i,j]
						non_opps = opportunity_matrices['1D'][i,j]
						
						snp_change_freqs[species_name][(sample_i, sample_j)].append((variant_type, freq_dict, syn_opps, non_opps, opp_dict)) #!
						'''
						# Calculate null version based on # of opportunities
						for prev_cohort in prev_cohorts:
							
							snv_freq_keys = snv_freq_map[prev_cohort].keys()
							
							L = snp_opportunity_matrix[i,j]
							L_snv = len(snv_freq_map[prev_cohort]) # A slight overestimate
							snv_fraction = L_snv*1.0/L
							num_bootstraps = 10
							
							for _ in range(num_bootstraps):
								
								if random() < snv_fraction: # Polymorphic site
									
									random_snv_location = snv_freq_keys[randint(0, len(snv_freq_keys))]
									f = 0 if random_snv_location in private_snv_map else snv_freq_map[prev_cohort][random_snv_location]
									rev_f = 1-f
									
									# Now add in probability weight
									snp_change_null_freqs[species_name][(sample_i, sample_j)].append((f, (1-f)*1.0/num_bootstraps))
									snp_change_null_freqs[species_name][(sample_i, sample_j)].append((1-f, f*1.0/num_bootstraps))
								
								else: # A truly invariant site
									
									snp_change_null_freqs[species_name][(sample_i, sample_j)].append((0, 1.0/num_bootstraps))
					
					for gene_change in gains:
						
						gene_name = gene_change[0]
						f = gene_freq_map[gene_name] if gene_name in gene_freq_map else 0
						gene_gain_freqs[species_name][(sample_i, sample_j)].append(f) #!
					
					for gene_change in losses:
						
						gene_name = gene_change[0]
						f = gene_freq_map[gene_name] if gene_name in gene_freq_map else 0
						gene_loss_freqs[species_name][(sample_i, sample_j)].append(f) #!
						num_bootstraps = 10
						fs = choice(gene_freq_values, size=num_bootstraps, p=gene_freq_weights)
						for f in fs:
							gene_loss_null_freqs[species_name][(sample_i, sample_j)].append(f)
				
				else: # Likely replacement and too many SNPs to store info for
					snp_changes[species_name][(sample_i, sample_j)] = num_snp_changes
					gene_changes[species_name][(sample_i, sample_j)] = (num_gains, num_losses)
			'''

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles" % ddir

# pickle.dump(snp_changes, open('%s/big_snp_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
# pickle.dump(gene_changes, open('%s/big_gene_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(snp_change_freqs, open('%s/snp_change_freqs_with_opps_%s.pkl' % (pdir, sweep_type), 'wb'))
# pickle.dump(snp_change_null_freqs, open('%s/snp_change_null_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
# pickle.dump(gene_gain_freqs, open('%s/gene_gain_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
# pickle.dump(gene_loss_freqs, open('%s/gene_loss_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
# pickle.dump(gene_loss_null_freqs, open('%s/gene_loss_null_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))

sys.stderr.write("Done!\n")
