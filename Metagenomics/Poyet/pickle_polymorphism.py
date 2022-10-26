# This script pickles polymorphism info

from utils import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import pylab, sys, numpy as np, random, math
from utils import temporal_changes_utils
from collections import defaultdict
import bz2

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
poyet_samples = su.get_sample_names('all')

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()
full_species_list = parse_midas_data.parse_species_list()

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
# Dictionary: sample(-species pair) -> within-sample polymorphism (average estimate)
sample_species_polymorphism_dict = defaultdict(dict)

# Same as above, but report (numerator, denominator) of polymorphism
sample_species_polymorphism_nd_dict = defaultdict(dict)
# ====================================================================

for species_name in full_species_list:
	
	print("Working on " + species_name)
	
	# Store within-host polymorphism for ALL samples
	_, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D']))
	
	for sample in poyet_samples:		
		
		if sample not in sfs_map: # TODO
			continue
		
		within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
		
		try:
			within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05) # 95% CI
		except:
			continue
		
		sample_species_polymorphism_dict[sample][species_name] = (within_rate_lower + within_rate_upper)/2.0
		
		sample_species_polymorphism_nd_dict[sample][species_name] = (within_sites, total_sites)

# Pickle things

import pickle
pickle_dir = "%s/pickles" % config.data_directory

pickle.dump(sample_species_polymorphism_dict, open("%s/sample_species_polymorphism_dict_full.pkl" % (pickle_dir), 'wb'))
pickle.dump(sample_species_polymorphism_nd_dict, open("%s/sample_species_polymorphism_nd_dict_full.pkl" % (pickle_dir), 'wb'))

# Store polymorphism and SNP change tables

analysis_dir = config.analysis_directory

with open("%s/Poyet_polymorphism_alpha_divs_full.csv" % analysis_dir, 'w') as outfile:
	fields = ','.join('sample_id subject_id species_name intfreq_sites_count total_sites_count polymorphism_rate Shannon_alpha_diversity'.split())
	outfile.write(fields + '\n')
	for sample in sample_species_polymorphism_dict:
		subject = sample_subject_map[sample]
		alpha_div = alpha_div_dict[sample]
		for species in sample_species_polymorphism_dict[sample]:
			wcount, tcount = sample_species_polymorphism_nd_dict[sample][species]
			prate = sample_species_polymorphism_dict[sample][species]
			fields_data = [sample, subject, species, wcount, tcount, prate, alpha_div]
			outfile.write(','.join([str(f) for f in fields_data]) + '\n')


