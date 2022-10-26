import config
import gzip, os
from collections import defaultdict

midas_dir = config.midas_directory

# ===========================================================================
# Loads set of genes in the MIDAS reference genome for a given species
# ===========================================================================

def load_reference_genes(species_name):
		
		rep_genome_filename = "%s/rep_genomes/%s/genome.features.gz" % (midas_dir, species_name)
		file = gzip.open(rep_genome_filename, 'r')		
		file.readline() # header
		
		reference_genes = []
		
		for line in file:
				items = line.split()
				gene_name = items[0].strip()
				reference_genes.append(gene_name)
		
		file.close()		
		
		return set(reference_genes)

# ===========================================================================
# Returns number of non-redundant genes in MIDAS pangenome for given species
# ===========================================================================

def get_reference_genome_size(species_name):
	return len(load_reference_genes(species_name))

# ===========================================================================
# Loads MIDAS pangenome map (genome -> gene -> centroids) for a given species
# ===========================================================================

def get_pangenome_map(species_name):
		
		gene_info_filename = '%s/pan_genomes/%s/gene_info.txt.gz' % (midas_dir, species_name)
		file = gzip.open(gene_info_filename, 'r')
		file.readline() # header
		
		pangenome_map = {}
		
		for line in file:
				items = line.split("\t")
				gene_id = items[0].strip()
				genome_id = items[1].strip()
				centroid_99 = items[2].strip()
				centroid_95 = items[3].strip()
				
				if genome_id not in pangenome_map:
						pangenome_map[genome_id] = {}
				
				pangenome_map[genome_id][gene_id] = (centroid_99, centroid_95)
		
		file.close()
		
		return pangenome_map

# ===========================================================================
# Loads set of non-redundant genes in MIDAS pangenome for given species
# where genes in same 99% identity cluster are considered redundant
# ===========================================================================

def load_pangenome_genes(species_name):
	
	pangenome_map = get_pangenome_map(species_name)
	
	non_redundant_genes = set()

	for genome in pangenome_map:
		for gene in pangenome_map[genome]:
			centroid_99 = pangenome_map[genome][gene][0]
			non_redundant_genes.add(centroid_99)
	
	return non_redundant_genes

# ===========================================================================
# Returns number of non-redundant genes in MIDAS pangenome for given species
# where genes in same 99% identity cluster are considered redundant
# ===========================================================================

def get_pangenome_size(species_name):
	return len(load_pangenome_genes(species_name))

# ===========================================================================
# Returns number of genomes in MIDAS pangenome for a given species
# ===========================================================================

def get_number_of_genomes(species_name):		
		return len(get_pangenome_map(species_name))

# ===========================================================================
# Returns list of MIDAS species (5926 species for midas_db_v1.2)
# ===========================================================================

def parse_species_list():
		
		species_directories = os.listdir("%s/pan_genomes" % midas_dir)
		
		species_names = []
		for potential_species_name in species_directories:
				if not potential_species_name.startswith('.'):
						species_names.append(potential_species_name)
		
		return species_names
		
# ===========================================================================
# The gene_ids in the pangenome list are the centroids of gene clusters.
# Sometimes the gene in the reference genome is not chosen as the centroid.
# This function creates a map between pangenome_centroids and genes in 
# reference genome (if it exists -- otherwise map to itself)
# ===========================================================================

def load_centroid_gene_map(species_name = None):
		
		if species_name == None:
				import parse_midas_data
				desired_species = parse_midas_data.parse_good_species_list()
		else:
				desired_species = [species_name]
		
		for species_name in desired_species:
				
				# First load reference genes
				reference_genes = load_reference_genes(species_name)
				
				# Next load pangenome centroids
				gene_info_filename = '%s/pan_genomes/%s/gene_info.txt.gz' % (midas_dir, species_name)
				gene_info_file = gzip.open(gene_info_filename, 'r')
				gene_info_file.readline() # header
		
				centroid_gene_map = {}
		
				for line in gene_info_file:
				
						items = line.split("\t") 
						gene_id = items[0].strip()
						centroid_id = items[3].strip() # 95% centroid
				
						if centroid_id not in centroid_gene_map:
								centroid_gene_map[centroid_id] = centroid_id
						
						if (gene_id in reference_genes) and (centroid_id not in reference_genes):
								centroid_gene_map[centroid_id] = gene_id
				
				gene_info_file.close()
		
		return centroid_gene_map

# ===========================================================================
# Returns set of genes in the MIDAS pangenome of a given species
# which are shared with other species (plus other shared genes which
# the given species may not have.)
# 
# The purpose is to ignore genes which have >= 95% sequence identity with
# at least one other gene in a different species' pangenome.
# ===========================================================================

def parse_midas_shared_genes(desired_species):
		
		midas_shared_genes = set()
		
		# Get mapping of pangenome centroids to reference genes
		centroid_gene_map = load_centroid_gene_map(desired_species)
		
		midas_db_shared_gene_filename = "%s/cross_species_centroids.txt.gz" % midas_dir
		file = gzip.open(midas_db_shared_gene_filename, "r")
		for line in file:
				items = line.split()
				
				# One name for a shared gene
				big_centroid = items[0].strip()
				midas_shared_genes.add(big_centroid)
				
				# Check if this species has any gene which is shared
				# If so, add it to the list of shared genes
				other_centroids = items[1].strip().split(",")
				for centroid in other_centroids:
						if centroid in centroid_gene_map:
								# Specifically, add the reference gene if possible
								midas_shared_genes.add(centroid_gene_map[centroid])
						
		return midas_shared_genes

# ===========================================================================
# Returns mapping of genes in the MIDAS pangenome of a given species
# to 'shared' genes in the pangenome (>= 95% sequence identity)
# ===========================================================================

def parse_midas_shared_genes_map():
		
		midas_shared_genes = defaultdict(set)
		
		midas_db_shared_gene_filename = "%s/cross_species_centroids.txt.gz" % midas_dir
		file = gzip.open(midas_db_shared_gene_filename, "r")
		for line in file:
				items = line.split()
				
				# One name for a shared gene
				big_centroid = items[0].strip()
				
				# Map it to set of similar genes
				other_centroids = items[1].strip().split(",")
				for centroid in other_centroids:
						midas_shared_genes[big_centroid].add(centroid)
		
		return midas_shared_genes

# ==========================================================================
# Find the representative genome id for a given species in MIDAS
# ==========================================================================

def representative_genome_id(desired_species_name):
	species_info = open("%s/species_info.txt" % midas_dir, 'r')
	species_info.readline() # header
	
	for line in species_info:
		species_id, rep_genome, count_genomes = line.strip().split("\t")
		if species_id == desired_species_name:
			return rep_genome

# ==========================================================================
# Get list of all genome ids corresponding to a species in MIDAS
# (i.e. every genome included in the pangenome)
# ==========================================================================

def get_ref_genome_ids(desired_species_name):
	
	genome_info = open("%s/genome_info.txt" % midas_dir)
	genome_info.readline() # header
	
	genome_ids = []
	for line in genome_info:
		genome_id, _, _, _, _, species_id = line.strip().split("\t")
		if desired_species_name == species_id:
			genome_ids.append(genome_id)
	
	return genome_ids

# ==========================================================================
# Return map: genome ID -> species_id			(for all species in MIDAS)
# ==========================================================================

def genome_ids_dictionary():
	
	genome_info = open("%s/genome_info.txt" % midas_dir)
	genome_info.readline() #header
	
	genome_dict = {}
	
	for line in genome_info:
		genome_id, _, rep_genome, _, _, species_id = line.strip().split("\t")
		genome_dict[genome_id]=[species_id, rep_genome]
	
	return genome_dict

# ==========================================================================
# Returns dictionary mapping all pangenome 99% centroids to their sequences
# ==========================================================================

def load_centroid_fasta(species_name):
		
		centroid_file = gzip.open("%s/pan_genomes/%s/centroids.ffn.gz" % (midas_dir, species_name), 'r')
		
		# Dictionary: centroid gene -> sequence
		centroid_fastas = defaultdict(str)
		
		for line in centroid_file:
			line = line.strip()
			if line[0] == '>':
				gene_name = line[1:len(line)]
			else:
				centroid_fastas[gene_name]+= line
		
		return centroid_fastas

# ==========================================================================
# METHODS FOR GETTING (95%) CLUSTER CENTROID INFO
# 	- load_centroid_gene_map	(maps centroid to reference gene if possible)
# 	- load_gene_centroid_map	(maps all pangenome genes to 95% centroids)
# 	- load_complete_centroid_gene_map	(maps centroid to all cluster genes)
# ==========================================================================

# ==========================================================================
# The gene_ids in the pangenome list are the centroids of gene clusters.
# Sometimes the gene in the reference genome is not chosen as the centroid.
# This function creates a map between pangenome_centroids and genes in 
# reference genome (if it exists)
# ==========================================================================

def load_centroid_gene_map(desired_species_name):
	
	# First load reference genes
	reference_genes = load_reference_genes(desired_species_name)
	
	# Load pangenome gene info from MIDAS db for this species
	gene_info_file = gzip.open("%s/pan_genomes/%s/gene_info.txt.gz" % (midas_dir, desired_species_name), 'r')	
	gene_info_file.readline() # remove header
	
	centroid_gene_map = {}
	
	for line in gene_info_file:
		gene_id, genome_id, c99, c95, c90, c85, c80, c75 = line.strip().split("\t") 
		centroid_id = c95 # choose 95% centroid
		
		if (gene_id in reference_genes) and (centroid_id not in reference_genes):
			centroid_gene_map[centroid_id] = gene_id
		elif centroid_id not in centroid_gene_map:
			centroid_gene_map[centroid_id] = centroid_id
	
	gene_info_file.close()	
	return centroid_gene_map

# ==========================================================================
# Returns map: pangenome gene name -> (reference corrected) centroid
# ==========================================================================

def load_gene_centroid_map(desired_species_name):
	
	# Maps centroid genes to corresponding reference genes if existing
	centroid_refgene_map = load_centroid_gene_map(desired_species_name)
	
	gene_info_file = gzip.open("%s/pan_genomes/%s/gene_info.txt.gz" % (midas_dir, desired_species_name), 'r')	
	gene_info_file.readline() # remove header
	
	gene_centroid_map = {}
	
	for line in gene_info_file:
		gene_id, genome_id, c99, c95, c90, c85, c80, c75 = line.strip().split("\t")
		ref_corrected_centroid = centroid_refgene_map[c95] # choose 95% centroid
		gene_centroid_map[gene_id] = ref_corrected_centroid
	
	gene_info_file.close()	
	return gene_centroid_map

# ==========================================================================
# Returns map: centroid -> list of ALL genes within the 95% cluster.
# ==========================================================================

def load_complete_centroid_gene_map(desired_species_name):
	
	gene_info_file = gzip.open("%s/pan_genomes/%s/gene_info.txt.gz" % (midas_dir, desired_species_name), 'r')	
	gene_info_file.readline() # remove header
	
	complete_centroid_gene_map = defaultdict(list)
	
	for line in gene_info_file:
		gene_id, genome_id, c99, c95, c90, c85, c80, c75 = line.strip().split("\t")
		complete_centroid_gene_map[c95].append(gene_id)
	
	gene_info_file.close()
	return complete_centroid_gene_map

# ==========================================================================
# Returns map: genus -> phylum
# ==========================================================================

def load_genus_phylum_map():
	
	genus_phylum_map = {}
	
	taxonomy_file = open("%s/genome_taxonomy.txt" % (midas_dir), 'r')
	taxonomy_file.readline() # remove header
	
	for line in taxonomy_file:
		genome_id, genome_name, taxon_id, kingdom, phylum, class_, order, family, genus, species, _, _ = line.strip('\n').split('\t')
		genus_phylum_map[genus] = phylum
	
	taxonomy_file.close()
	return genus_phylum_map

# ==========================================================================
# Returns map: genus / family / order -> phylum
# ==========================================================================

def load_gfo_phylum_map():
	
	genus_phylum_map = {}
	
	taxonomy_file = open("%s/genome_taxonomy.txt" % (midas_dir), 'r')
	taxonomy_file.readline() # remove header
	
	for line in taxonomy_file:
		genome_id, genome_name, taxon_id, kingdom, phylum, class_, order, family, genus, species, _, _ = line.strip('\n').split('\t')
		# Manually correct for Bilophila
		if genus == 'Bilophila':
			genus_phylum_map[genus] = 'Proteobacteria'
		else:
			genus_phylum_map[genus] = phylum
		genus_phylum_map[family] = phylum
		genus_phylum_map[order] = phylum
	
	taxonomy_file.close()
	return genus_phylum_map

# ===========================================================================
# Outputs information on number of genomes per species, why not?
# ===========================================================================

if __name__=='__main__':
		
		import parse_midas_data
		good_species_list = parse_midas_data.parse_good_species_list()
		for species_name in good_species_list:
				num_genomes = get_number_of_genomes(species_name)
				print("%s: %i" % (species_name, num_genomes))


