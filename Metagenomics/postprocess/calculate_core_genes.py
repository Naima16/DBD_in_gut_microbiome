import numpy, sys, gzip, os.path, os
from utils import config, midas_db_utils, parse_midas_data, core_gene_utils

core_genes_directory = "%s/core_genes/" % config.data_directory
external_core_genes_directory = "%s/core_genes/external/" % config.data_directory

default_external_shared_gene_filename					= "%s/shared_genes.txt.gz" % external_core_genes_directory
default_external_core_gene_filename						= "%s/core_genes.txt.gz" % external_core_genes_directory
default_external_stringent_core_gene_filename = "%s/core_genes_stringent.txt.gz" % external_core_genes_directory
default_external_gene_freq_template						= (external_core_genes_directory + "%s_gene_freqs.txt.gz")

default_shared_gene_filename					= "%s/shared_genes.txt.gz" % core_genes_directory
default_core_gene_filename						= "%s/core_genes.txt.gz" % core_genes_directory
default_stringent_core_gene_filename	= "%s/core_genes_stringent.txt.gz" % core_genes_directory
default_gene_freq_template						= (core_genes_directory+"%s_gene_freqs.txt.gz")

# ===========================================================================
# Computes lists of shared genes, core genes, and "stringent" core genes for
# each species, as well as prevalences of each gene
# ===========================================================================

os.system('mkdir -p %s' % core_genes_directory)
os.system('mkdir -p %s' % external_core_genes_directory)

pangenome_species = parse_midas_data.parse_good_species_list()

cmin = config.core_genome_min_copynum # 0.3
cmax = config.core_genome_max_copynum	# 3
shared_cmin = config.shared_genome_min_copynum # 3

min_good_fraction = config.core_genome_min_prevalence # 0.9
min_coverage = 5 # (for assessing core genome, we'll use a lower coverage value than when we look at real changes)

output_file = gzip.GzipFile(default_core_gene_filename,"w")
stringent_output_file = gzip.GzipFile(default_stringent_core_gene_filename,"w")
shared_output_file = gzip.GzipFile(default_shared_gene_filename,"w")

for species_name in pangenome_species:
		
		# Load reference genes
		sys.stderr.write("Loading genes on reference genome for %s..\n" % species_name)
		reference_genes = midas_db_utils.load_reference_genes(species_name)
		sys.stderr.write("Done!\n")
		
		# Load shared genes
		sys.stderr.write("Loading shared genes from midas db for %s..\n" % species_name)
		midas_shared_genes = midas_db_utils.parse_midas_shared_genes(species_name)
		sys.stderr.write("Done!\n")
		
		# Load gene coverage information for species_name
		sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
		gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
		sys.stderr.write("Done!\n")	 
		
		# This turns to true if not enough high coverage samples
		bad_pangenome_data = False
		
		if len(marker_coverages)==0: # TODO: why would this happen?
				bad_pangenome_data = True
		else:				 
				
				high_coverage_idxs = (marker_coverages>=min_coverage)

				if high_coverage_idxs.sum() < 0.5:
						bad_pangenome_data = True
		
		if bad_pangenome_data:
				
				# Just use reference genes
				sys.stderr.write("Bad pangenome data for %s!\n" % species_name)
				shared_gene_names = sorted(midas_shared_genes)
				core_gene_names = sorted(reference_genes - midas_shared_genes)
				stringent_gene_names = sorted(reference_genes - midas_shared_genes)
	
		else:		 
		
				gene_names = numpy.array(gene_names)
				
				# Only look at samples with marker coverage >5
				
				gene_samples = gene_samples[high_coverage_idxs]
				marker_coverages = marker_coverages[high_coverage_idxs]
				gene_depth_matrix = gene_depth_matrix[:,high_coverage_idxs]
				
				# Note that samples with marker coverage of 0 end up
				# with dummy marker coverage of 1, so that copy number
				# is equal to depth for all genes.
				
				# These cases are weird because we would expect 0 depth
				# for all genes of, say, B. vulgatus in a sample,
				# if its marker gene coverage is 0 (suggesting that
				# the species is not present at all).
				# So if we get nonzero depth for those other genes,
				# we leave it be -- but probably bioinformatics error.
				
				gene_copynum_matrix = gene_depth_matrix.astype('float')/(marker_coverages+(marker_coverages==0))
				
				# Further restrict to samples with reasonable copy numbers for most genes
				
				good_sample_idxs = core_gene_utils.get_good_pangenome_samples(marker_coverages, gene_copynum_matrix, species_name)
				bad_sample_idxs = numpy.logical_not(good_sample_idxs)
				sys.stderr.write("%d bad samples!\n" % bad_sample_idxs.sum())
				
				gene_samples = gene_samples[good_sample_idxs]
				marker_coverages = marker_coverages[good_sample_idxs]
				gene_copynum_matrix = gene_copynum_matrix[:,good_sample_idxs]
				
				# Find out which genes are in the reference genome, and which are shared
				
				reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
				midas_shared_idxs = numpy.array([gene_name in midas_shared_genes for gene_name in gene_names])
				
				# These are genes that have coverage >=3x normal in some sample. This are candidates for being linked to another species.
				# (they could also be multi-copy genes, but we can't look at much on these genes anyway, so might as well toss them out)
				
				metagenome_shared_idxs = ((gene_copynum_matrix>shared_cmin).sum(axis=1)>0.5)
				
				# Now union with those shared genes we identified from midas db
				
				shared_idxs = numpy.logical_or(metagenome_shared_idxs, midas_shared_idxs)
				non_shared_idxs = numpy.logical_not(shared_idxs)
				
				shared_gene_names = gene_names[shared_idxs]
				
				# Calculating good genes.
				# Note that * is logical AND, so we are getting, for each gene,
				# the number of samples that the gene has intermediate copynum in,
				# dividing by the total number of samples, and checking if this fraction
				# is greater than 90%.
				
				good_idxs = (((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)*1.0/len(marker_coverages) >= min_good_fraction)
				
				# Define "core genes" as those which are in the reference genome,
				# not shared, and have reasonably copynum in the majority of samples
				# (samples themselves have been filtered for high coverage and reasonable
				# copy numbers for most of its genes)
				
				core_gene_idxs = good_idxs*reference_gene_idxs*non_shared_idxs
				core_gene_names = gene_names[core_gene_idxs]
				
				# Summary statistics on how many genes are shared (midas or metagenome) or not
				
				num_metagenome_and_midas = numpy.logical_and(midas_shared_idxs, metagenome_shared_idxs).sum()
				num_metagenome_only = numpy.logical_and(metagenome_shared_idxs, numpy.logical_not(midas_shared_idxs)).sum()
				num_midas_only = numpy.logical_and(midas_shared_idxs, numpy.logical_not(metagenome_shared_idxs)).sum()
				num_metagenome_or_midas = shared_idxs.sum()
				num_remaining = non_shared_idxs.sum()
				num_reference_remaining = (non_shared_idxs*reference_gene_idxs).sum()
				num_core = core_gene_idxs.sum()
				
				print "%s %d %d %d %d %d %d %d" % (species_name, num_metagenome_and_midas, num_metagenome_only, num_midas_only, num_metagenome_or_midas, num_remaining, num_reference_remaining, num_core)
				
				# Measure frequencies and output them
				
				# For each gene, number of samples that the gene has intermediate copynum in
				gene_prevalence_numerators = ((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)
				
				# For each gene, number of samples that gene has not-too-high copynum in
				gene_prevalence_denominators = ((gene_copynum_matrix<=cmax).sum(axis=1))
				
				# Genes that are not shared + have intermediate copynum in at least one sample
				good_prevalence_idxs = (gene_prevalence_numerators>0.5)*(gene_prevalence_denominators>0.5)*non_shared_idxs
				gene_prevalence_names = gene_names[good_prevalence_idxs]
				
				# Compute prevalence of such "good" genes
				gene_prevalences = gene_prevalence_numerators[good_prevalence_idxs]*1.0/gene_prevalence_denominators[good_prevalence_idxs]
				
				# Save "good" gene prevalence information
				gene_freq_output_file = gzip.GzipFile(default_gene_freq_template % species_name,"w")
				for gene_name, f in zip(gene_prevalence_names, gene_prevalences):
						gene_freq_output_file.write("%s %g\n" % (gene_name, f))
				gene_freq_output_file.close()

				# Calculating good genes w/ stringent definition (100%)
				# That is, for each gene, if it has < 0.05 copynum in at least one sample,
				# it is considered not good.
				
				bad_idxs = (gene_copynum_matrix<config.gainloss_max_absent_copynum).sum(axis=1) > 0.5
				good_idxs = numpy.logical_not(bad_idxs)
				stringent_gene_names = gene_names[good_idxs*reference_gene_idxs*non_shared_idxs]
				
				sys.stderr.write("%d stringent core genes out of %d\n" % (len(stringent_gene_names), len(gene_names)))		 
		
		# Write output to file!
		shared_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in shared_gene_names])))
		output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in core_gene_names])))
		stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in stringent_gene_names])))

shared_output_file.close()
output_file.close()
stringent_output_file.close()
