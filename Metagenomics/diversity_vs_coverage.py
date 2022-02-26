from utils import sample_utils as su, parse_midas_data, config
import bz2, sys
import numpy as np
import math

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Samples
hmp_samples = su.get_sample_names('all')

# Species list
species_list = parse_midas_data.parse_species_list()

# Relative abundance file
relab_fpath = "%s/species/relative_abundance.txt.bz2" % (config.data_directory)
relab_file = open(relab_fpath, 'r')
decompressor = bz2.BZ2Decompressor()
raw = decompressor.decompress(relab_file.read())
data = [row.split('\t') for row in raw.split('\n')]
data.pop() # Get rid of extra element due to terminal newline
header = data[0]

# Generate alpha diversity dictionary
alpha_div_dict = {}
for i in range(1, len(header)):
	acc = 0
	for row in data[1:]:
		rel_ab = float(row[i])
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
	alpha_div_dict[header[i]] = (acc*-1)

# Get average marker gene coverage for each sample-species pair
cov_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()

# Sum up average marker gene coverage across all species
agg_cov_matrix = np.zeros(len(samples))

# Average marker gene coverages for each species
species_cov_dict = {}

for species_covs, cur_species in zip(cov_matrix, species):
	agg_cov_matrix += species_covs
	species_cov_dict[cur_species] = species_covs

ordered_alpha_divs = [alpha_div_dict[s] for s in samples]

from matplotlib import pyplot as plt

fig, ax = plt.subplots()
ax.plot(ordered_alpha_divs, agg_cov_matrix, '.')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Average coverage of marker genes summed over all species")
ax.set_title("Sample coverage vs. species diversity")
fig.savefig('%s/coverage_vs_diversity.png' % (config.analysis_directory),bbox_inches='tight')

# Plot 2
cov_sorted_species_covs = sorted(species_cov_dict.values(), key=lambda covs: covs.mean())

all_species = species_cov_dict.keys()
all_species_covs = [species_cov_dict[s] for s in all_species]

fig, ax = plt.subplots()

ax.boxplot(cov_sorted_species_covs)
fig.savefig('%s/coverage_by_species_vs_diversity.png' % (config.analysis_directory),bbox_inches='tight')