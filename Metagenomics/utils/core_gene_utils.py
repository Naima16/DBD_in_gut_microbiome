import numpy, os.path, gzip
import config, parse_midas_data

# ===========================================================================
# This module (core_gene_utils) contains the following utilities:
# 
# 	parse_core_genes
# 	parse_shared_genes
# 	parse_non_shared_reference_genes
# 	parse_non_shared_pangenome_genes
# 	get_good_pangenome_samples
# 	parse_gene_freqs
#
# ===========================================================================

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
# Returns set of core genes for specified set of species
# ===========================================================================

def parse_core_genes(desired_species_name = None, core_gene_filename=default_core_gene_filename, external_core_gene_filename=default_external_core_gene_filename, external_filtering=True):
		
		core_genes = set() # Core genes for the specified species
		core_gene_file = gzip.GzipFile(core_gene_filename,"r")
		for line in core_gene_file:
				
				items = line.split(":")
				if len(items)<2:
						continue
						
				species_name = items[0].strip()
				gene_names = [subitem.strip() for subitem in items[1].split(",")]
				
				if (species_name == desired_species_name) or (desired_species_name == None):
						core_genes.update(gene_names)
						
		core_gene_file.close() 
		
		external_core_genes = set()
		if os.path.isfile(external_core_gene_filename):
				
				external_core_gene_file = gzip.GzipFile(external_core_gene_filename,"r")
				
				for line in external_core_gene_file:
				
						items = line.split(":")
						if len(items)<2:
								continue
						
						species_name = items[0].strip()
						gene_names = [subitem.strip() for subitem in items[1].split(",")]
				
						if (species_name == desired_species_name) or (desired_species_name == None):
								external_core_genes.update(gene_names)
						
				external_core_gene_file.close() 
		
		if external_filtering and len(external_core_genes)>0:
				# some externally provided core genes
				core_genes = (core_genes & external_core_genes)
				
		return core_genes

# ===========================================================================
# Returns set of shared genes for specified set of species
# ===========================================================================

def parse_shared_genes(desired_species_name = None, shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering = True):
		
		shared_genes = set()
		shared_gene_file = gzip.GzipFile(shared_gene_filename,"r")
		for line in shared_gene_file:
				
				items = line.split(":")
				if len(items)<2:
						continue
						
				species_name = items[0].strip()
				gene_names_str = items[1].strip()
				if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
						gene_names = []
				else:
						gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
				
				if (species_name==desired_species_name) or (desired_species_name==""):
						shared_genes.update(gene_names)
						
		shared_gene_file.close() 
		
		external_shared_genes = set()
		if os.path.isfile(external_shared_gene_filename):
				
				external_shared_gene_file = gzip.GzipFile(external_shared_gene_filename,"r")
						
				for line in external_shared_gene_file:
				
						items = line.split(":")
						if len(items)<2:
								continue
						
						species_name = items[0].strip()
						gene_names_str = items[1].strip()
						if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
								gene_names = []
						else:
								gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
				
						if (species_name==desired_species_name) or (desired_species_name==""):
								external_shared_genes.update(gene_names)
						
				external_shared_gene_file.close() 
		
		if external_filtering and len(external_shared_genes)>0:
				# some externally provided core genes
				shared_genes = (shared_genes | external_shared_genes)
				
		return shared_genes

# ===========================================================================
# Returns set of reference genes which are not shared
# ===========================================================================

def parse_non_shared_reference_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
		
		from utils import parse_midas_data
		
		shared_genes = parse_shared_genes(desired_species_name, shared_gene_filename, external_shared_gene_filename, external_filtering)
		reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
		non_shared_reference_genes = set(reference_genes)-shared_genes
		
		return non_shared_reference_genes

# ===========================================================================
# Returns set of pangenome genes which are not shared
# ===========================================================================

def parse_non_shared_pangenome_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
		
		from utils import parse_midas_data
		
		shared_genes = parse_shared_genes(desired_species_name, shared_gene_filename, external_shared_gene_filename, external_filtering)
		pangenome_genes, pangenome_centroid_genes = parse_midas_data.load_pangenome_genes(desired_species_name)
		# TODO: Not sure if I should be using the first or second
		non_shared_pangenome_genes = set(pangenome_centroid_genes)-shared_genes
		return non_shared_pangenome_genes

# ===========================================================================
# Returns indices for samples which have enough present genes (copy number
# exceeds a low threshold) of which not too many are high-copynum
# ===========================================================================

def get_good_pangenome_samples(marker_coverages, gene_copynum_matrix, species_name):
		
		cmin = config.core_genome_min_copynum
		cmax = config.core_genome_max_copynum	 
		
		# For each sample, get number of genes which are "present" (copynum > 0.3)
		# and which are "high" (copynum > 3)
		num_present_genes = (gene_copynum_matrix>cmin).sum(axis=0)
		num_high_genes = (gene_copynum_matrix>cmax).sum(axis=0)
		
		# Get number of reference genes
		# TODO: does not match number of genes in copynum matrix
		num_reference_genes = len(parse_midas_data.load_reference_genes(species_name))
		
		# Want at least 30% of all reference genes to be present, and want no more
		# than 30% of present genes to be high
		min_present_genes = 0.3*num_reference_genes
		max_high_genes = 0.3*num_present_genes
		
		good_sample_idxs = (num_present_genes>min_present_genes)*(num_high_genes<max_high_genes)
		
		return good_sample_idxs

# ===========================================================================
# Returns gene frequency map (gene name -> prevalence)
# ===========================================================================

def parse_gene_freqs(desired_species_name, use_external=False):
		
		if use_external:
				filename_template = default_external_gene_freq_template
		else:
				filename_template = default_gene_freq_template
		
		
		filename = filename_template % (desired_species_name)
		if not os.path.isfile(filename):
				return None
				
		file = gzip.open(filename,"r")
		gene_freq_map = {}
		for line in file:
				items = line.split()
				gene_name = items[0]
				f = float(items[1])
				gene_freq_map[gene_name] = f
		file.close()
		
		return gene_freq_map