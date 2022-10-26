# from parse_midas_data

# ===========================================================================
# Loads marker gene coverages (sample -> marker gene coverage) for a species
# ===========================================================================

def load_marker_gene_coverages(species_name):
	
	# genes_summary.txt: use to get marker gene coverage
	# i.e. median read-depth across 15 universal single copy genes
	gene_summary_file = file("%s/genes/%s/genes_summary.txt" % (data_dir, species_name),"r")
	
	marker_coverage_map = {}
	gene_summary_file.readline() # header
	marker_coverage_samples = []
	marker_coverages = []
	for summary_line in gene_summary_file:
			items = summary_line.split()
			sample = items[0].strip()
			marker_coverage = float(items[5])
			marker_coverage_samples.append(sample)
			marker_coverages.append(marker_coverage)
	
	gene_summary_file.close()
	
	marker_coverage_samples = su.parse_merged_sample_names(marker_coverage_samples)
	marker_coverage_map = {s: cov for s, cov in zip(marker_coverage_samples, marker_coverages)}
	
	return marker_coverage_map

# ===========================================================================
# Loads metaphlan2 genes (that are present in the reference genome)
# returns a list of metaphlan2 genes
# ===========================================================================

def load_metaphlan2_genes(desired_species_name):
	gene_file = open("%s/metaphlan2_genes/%s_metaphlan2_genes_mapped.txt" % (data_dir, desired_species_name), 'r')
	
	reference_genes = load_reference_genes(desired_species_name)
	
	metaphlan2_genes=[]
	for line in gene_file:
			gene_name = line.strip()
			
			if gene_name in reference_genes:
					metaphlan2_genes.append(gene_name)
			else:
					pass
					#print gene_name

	gene_file.close()		 
	
	return set(metaphlan2_genes)

def calculate_relative_depth_threshold_map(sample_coverage_histograms, samples, min_nonzero_median_coverage=5, lower_factor=0.5, upper_factor=2):
	
	# returns map of sample name: coverage threshold
	# essentially filtering out samples whose marker depth coverage
	# does not exceed the average coverage threshold
	
	depth_threshold_map = {}
	for i in xrange(0,len(samples)):
			
			# Check if coverage distribution meets certain requirements
			is_bad_coverage_distribution = False
			
			# First check if passes median coverage requirement
			nonzero_median_coverage = stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histograms[i])
			if round(nonzero_median_coverage) < min_nonzero_median_coverage:
					is_bad_coverage_distribution=True
	
			# Passed median coverage requirement
			# Now check whether a significant number of sites fall between lower and upper factor. 
			lower_depth_threshold = floor(nonzero_median_coverage*lower_factor)-0.5 # why is 0.5 being added/subtracted? NRG
			upper_depth_threshold = ceil(nonzero_median_coverage*upper_factor)+0.5
	
			depths, depth_CDF = stats_utils.calculate_CDF_from_histogram(sample_coverage_histograms[i])
			# remove zeros
			if depths[0]<0.5:
					depth_CDF -= depth_CDF[0]
					depth_CDF /= depth_CDF[-1]
			
			fraction_in_good_range = depth_CDF[(depths>lower_depth_threshold)*(depths<upper_depth_threshold)].sum()
	
			if fraction_in_good_range < 0.6: #where does 0.6 come from? NRG
					is_bad_coverage_distribution=True
					
			if is_bad_coverage_distribution:
					lower_depth_threshold = 1000000001
					upper_depth_threshold = 1000000001
			
			depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
			
	return depth_threshold_map


def calculate_absolute_depth_threshold_map(species_coverage_vector, samples, avg_depth_threshold=20, site_depth_threshold=15):
		
		# returns map of sample name: coverage threshold
		# essentially filtering out samples whose marker depth coverage
		# does not exceed the average coverage threshold
		
		depth_threshold_map = {}
		for i in xrange(0,len(samples)):
				
				if species_coverage_vector[i]<avg_depth_threshold:		
						lower_depth_threshold=1000000001
				else:
						lower_depth_threshold=site_depth_threshold
		
				upper_depth_threshold = 1000000001
				depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
				
		return depth_threshold_map