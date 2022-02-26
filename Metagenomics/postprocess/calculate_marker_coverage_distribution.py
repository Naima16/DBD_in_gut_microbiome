from utils import parse_midas_data, config
import sys, bz2
import numpy

if len(sys.argv) > 1:
	species_name = sys.argv[1]
else:
	sys.stderr.write("Usage: marker_coverage_distribution.py species")

sys.stderr.write("Calculating marker gene coverage distribution for %s...\n" % species_name)

# ==========================================================================
# Creates a coverage_distribution.txt.bz2 file in the species snp directory
# 
# In this file, rows are samples, columns are D,count pairs
# samples are guaranteed to be same order as snps_depth.txt.bz2 file
#
# Also creates a gene_coverage.txt.bz2 file in the species snp directory
# In this file, rows are genes, columns are samples, 
# entries are avg coverage of that gene for that sample 
# ==========================================================================

# Specify variants: use all types for most information
allowed_variant_types = set(["1D","2D","3D","4D"])

# Load list of marker genes
marker_genes = parse_midas_data.load_marker_genes(species_name)

# Read depth information straight from snps_depth.txt. No preprocessing
depth_file = bz2.BZ2File("%s/snps/%s/snps_depth.txt.bz2" % (config.data_directory, species_name), 'r')
depth_header = depth_file.readline() # header

# To get information about which gene it is in
info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (config.data_directory, species_name), 'r')
info_header = info_file.readline() # header

# Header of snps_depth: all samples that have this species
samples = depth_header.strip().split()[1:]

# ==============================================================================
# Dictionary: marker gene -> coverages of all sites with known variant type
# note that values are lists of coverages across all samples
marker_coverages = {gene_name: [] for gene_name in marker_genes}

# Dictionary: marker gene -> location where it starts
marker_gene_starts = {gene_name: -1 for gene_name in marker_genes}

# Dictionary: marker gene -> locations of all sites with known variant type
# where location is relative to the start location of gene w.r.t. whole genome
marker_locations = {gene_name: [] for gene_name in marker_genes}
# ==============================================================================

num_sites_processed = 0
total_marker_sites = 0

while True:
	
	num_sites_processed += 1		
	if num_sites_processed % 100000 == 0:
		sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed/1000))
	
	# load next lines
	depth_line = depth_file.readline()
	info_line = info_file.readline()
	
	# quit if file has ended
	if depth_line == "":
		break
	
	# parse site info
	site, _, _, _, _, variant_type, gene_name, _, _ = info_line.split('\t')
	
	if (variant_type not in allowed_variant_types) or (gene_name not in marker_genes):
		continue
	
	items = depth_line.strip().split('\t')
	location = long(items[0].split("|")[1])
	depths = numpy.array([long(item) for item in items[1:]])
	
	if marker_gene_starts[gene_name] < 0:
		marker_gene_starts[gene_name] = location
	
	marker_locations[gene_name].append(location - marker_gene_starts[gene_name])
	marker_coverages[gene_name].append(depths)

depth_file.close()
info_file.close()

# Now write output!

# First write genome-wide coverage distribution
output_fname = "%s/snps/%s/marker_coverage_distribution.txt.bz2" % (config.data_directory, species_name)
output_file = bz2.BZ2File(output_fname, 'w')

output_file.write("SampleID,GeneID\tD,n(D) ...\n")

# For each sample that has SNPS output for this species,
# 	for each marker gene,
# 		write list of (loc, cov) where loc is location within the gene

for sample_idx in range(len(samples)):
	sample = samples[sample_idx]
	for gene_name in marker_genes:
		output_file.write("%s,%s" % (sample,gene_name))
		
		for i in range(len(marker_locations[gene_name])):
			loc = marker_locations[gene_name][i]
			depth = marker_coverages[gene_name][i][sample_idx]			
			output_file.write('\t')
			output_file.write("%d,%d" % (loc, depth))
		
		output_file.write("\n")

output_file.close()
