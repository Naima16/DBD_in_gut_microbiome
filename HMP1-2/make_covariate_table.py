from utils import config, sample_utils as su
from collections import defaultdict
import math, bz2

data_dir = config.data_directory

# ========================================
# Make a table with the following fields:
# ========================================
# Sample ID
# Alpha diversity
# Richness (number of species)
# Number of reads (after MIDAS species)
# Original number of reads (from FASTQ)
# Total average coverage of marker genes
# Total median coverage of marker genes
# ========================================

# Relative abundance file (for alpha diversity, richness)
relab_fpath = "%s/species/relative_abundance.txt.bz2" % data_dir
relab_file = bz2.BZ2File(relab_fpath, 'r')
data = [row.strip().split('\t') for row in relab_file]
header = data[0]

alpha_div_dict = {}
richness_dict = {}

for i in range(1, len(header)):
	acc = 0
	richness_acc = 0
	
	for row in data[1:]:
		rel_ab = float(row[i])
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
			richness_acc += 1
	
	alpha_div_dict[header[i]] = (acc*-1)
	richness_dict[header[i]] = richness_acc

# Original read counts file (for number of reads in FASTQ files)
orig_read_count_dict = su.parse_sample_read_count_map()
'''
samples = orig_read_count_dict.keys()
for sample in samples:
	if orig_read_count_dict[sample] > 20000000:
		orig_read_count_dict[sample] = 20000000
'''
# Read counts file (for number of reads mapped to some species)
readc_fpath = "%s/species/count_reads.txt.bz2" % data_dir
readc_file = bz2.BZ2File(readc_fpath, 'r')
data = [row.strip().split('\t') for row in readc_file]
header = data[0]

read_count_dict = defaultdict(int)

for row in data[1:]: # Loops over species
	for i in range(1, len(header)):
		read_count = int(row[i])
		read_count_dict[header[i]] += read_count

# Coverage file: average read-depth of 15 marker genes per species
# (total bp of mapped reads/total bp of 15 marker-genes)
cov_file = bz2.BZ2File("%s/species/coverage.txt.bz2" % data_dir, 'r')
samples = cov_file.readline().strip().split()[1:] # header

total_avg_marker_cov_dict = defaultdict(float)

for line in cov_file:
	covs = line.strip().split()[1:]
	for i in range(len(samples)):
		total_avg_marker_cov_dict[samples[i]] += float(covs[i])

# Finally, get median coverage of marker gene summed over species
total_med_marker_cov_dict = defaultdict(float)

# Get list of species that MIDAS analyzed genes for...
# (There are 194)
with open("%s/genes/species_genes.txt" % data_dir, 'r') as file:
	species_names = [line.strip() for line in file]

# Genes summary file: contains median read-depth across 15 marker genes
# Note that median marker coverage is used to estimate copy number
for species in species_names:
	gsumm_fpath = "%s/genes/%s/genes_summary.txt" % (data_dir, species)
	gsumm_file = open(gsumm_fpath, 'r')
	gsumm_file.readline() # remove header
	for line in gsumm_file:
		sample, _, _, _, _, med_cov = line.strip().split()
		total_med_marker_cov_dict[sample] += float(med_cov)

# Sample subject map
sample_subject_map = su.parse_sample_subject_map()

# Write to sample_covariates.txt (under analysis directory)
output_file = open("%s/sample_covariates.txt" % config.analysis_directory, 'w')
fields = ["sample_id", "subject_id", "shannon_diversity", "species_richness", "total_reads_MIDAS", "total_reads_orig", "total_avg_marker_coverage", "total_med_marker_coverage"]
output_file.write('\t'.join(fields) + '\n')

samples = su.get_sample_names(remove_c = False)

for sample in samples:	
	
	sample_no_c = sample[:-1] if sample[-1] == 'c' else sample
	
	field_data = [sample_no_c, sample_subject_map[sample_no_c], alpha_div_dict[sample], richness_dict[sample], read_count_dict[sample], orig_read_count_dict[sample], total_avg_marker_cov_dict[sample], total_med_marker_cov_dict[sample]]
	output_file.write('\t'.join([str(f) for f in field_data]) + '\n')

output_file.close()
