import sys, bz2
import numpy
from collections import defaultdict
from utils import config

if len(sys.argv) > 1:
	species_name = sys.argv[1]
else:
	sys.stderr.write("Usage: coverage_distribution.py species")

sys.stderr.write("Calculating coverage distribution for %s...\n" % species_name)

# ===========================================================================
# Creates coverage_distribution.txt.bz2 file in species snps directory
# 
# In this file, rows are samples, columns are D,count pairs
# samples are guaranteed to be same order as snps_depth.txt.bz2 file
#
# Also creates a gene_coverage.txt.bz2 file in the species snp directory
# In this file, rows are genes, columns are samples, 
# entries are avg coverage of that gene for that sample 
# ===========================================================================

# These are the default MIDAS parameters. Used to ensure consistency.
# Only sites that pass the prevalence threshold in terms of having at least
# the prevalence_min_coverage can pass, stored in sample_depth_histograms

prevalence_threshold = 0.95
prevalence_min_coverage = 3

# Use all types of sites to include most information
allowed_variant_types = set(["1D","2D","3D","4D"]) 

snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)

# snps_depth.txt: number of reads mapped to genomic site per sample
depth_file = bz2.BZ2File("%s/snps_depth.txt.bz2" % snps_dir, 'r')

# snps_info.txt: metadata for genomic site: mean freq, mean depth, allele props, etc.
info_file = bz2.BZ2File("%s/snps_info.txt.bz2" % snps_dir, 'r')

# Get samples from depth file		
samples = depth_file.readline().split()[1:]

info_file.readline() # remove header

# sample_depth_histograms: (sample -> site -> depth)
sample_depth_histograms = {sample: defaultdict(int) for sample in samples} # stores only prevalent sites
full_sample_depth_histograms = {sample: defaultdict(int) for sample in samples} # stores all sites

gene_total_depths = {} # gene -> total depth across all considered sites
gene_total_sites = {} # gene -> number of sites (with known variant type)

num_sites_processed = 0

while True:
		
		# load next lines
		depth_line = depth_file.readline()
		info_line = info_file.readline()
		
		# quit if file has ended
		if depth_line=="":
				break
		
		# parse site info
		info_items = info_line.split('\t')
		variant_type = info_items[5]
				
		# make sure it is a site with known variant type
		if variant_type not in allowed_variant_types:
				continue
		
		gene_name = info_items[6]
		
		depth_items = depth_line.split()
		depths = numpy.array([long(item) for item in depth_items[1:]])
		
		# Manual prevalence filter
		if (depths>=prevalence_min_coverage).sum()*1.0/len(depths) >= prevalence_threshold: 
				# Add to genome-wide depth distribution
				for sample, D in zip(samples,depths):
						sample_depth_histograms[sample][D] += 1
		
		for sample, D in zip(samples,depths):
				full_sample_depth_histograms[sample][D] += 1
		
		# Add to gene-specific avg depth (not subject to prevalence filter)
		if gene_name not in gene_total_depths:
				gene_total_depths[gene_name] = numpy.zeros(len(samples))*1.0
				gene_total_sites[gene_name] = 0.0
				
		gene_total_depths[gene_name] += depths
		gene_total_sites[gene_name] += 1
		
		num_sites_processed += 1
		
		if num_sites_processed%100000==0:
				sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed/1000))

depth_file.close()

# Now write output!

# ===========================================================================
# First write (prevalence-filtered) genome-wide coverage distribution
# 
# Details:
# 	- rows are samples (same as snps_depth.txt samples)
# 	- columns are depths (varies by sample)
# 	- items are number of sites which have that depth, genome-wide
# 
# N.B. we are only considering sites with known variant type, not
# all sites in the genome. Further filtered for having enough coverage
# in at least 95% of all samples.
# ===========================================================================

output_file = bz2.BZ2File("%s/coverage_distribution.txt.bz2" % snps_dir, 'w')
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
		output_file.write("\n")
		output_file.write("\t".join([sample]+["%d,%d" % (D,sample_depth_histograms[sample][D]) for D in sorted(sample_depth_histograms[sample].keys())]))

output_file.close()

# ===========================================================================
# Write unfiltered genome-wide coverage distribution
# 
# Same form as above, except consider all sites with known variant type.
# ===========================================================================

output_file = bz2.BZ2File("%s/full_coverage_distribution.txt.bz2" % snps_dir, 'w')
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
		output_file.write("\n")
		output_file.write("\t".join([sample]+["%d,%d" % (D, full_sample_depth_histograms[sample][D]) for D in sorted(full_sample_depth_histograms[sample].keys())]))

output_file.close()

# ===========================================================================
# Then write gene-specific coverages
# 
# Details:
# 	- rows are genes
# 	- columns are samples
# 	- items are average depth of all known sites in that gene
# 
# Note: gene_total_sites[gene_name] == 0 should never happen...
# so I got rid of it...
# ===========================================================================

output_file = bz2.BZ2File("%s/gene_coverage.txt.bz2" % snps_dir, 'w')
output_file.write("\t".join(["Gene"]+samples)) # Header line
for gene_name in sorted(gene_total_depths.keys()):
		avg_depths = gene_total_depths[gene_name]/(gene_total_sites[gene_name])
		output_file.write("\n")
		output_file.write("\t".join([gene_name]+["%0.1f" % D for D in avg_depths]))

output_file.close()

# Done!