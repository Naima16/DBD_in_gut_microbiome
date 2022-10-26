from collections import defaultdict
from utils import sample_utils as su
import bz2

'''
sample = '700175181'
gene_fam_fpath = '%s_1_pathabundance.tsv' % (sample)

pathway_data = {}
pathway_species_data = defaultdict()

with open(gene_fam_fpath, 'r') as file:
	header = file.readline()
	for line in file:
		gene_fam, rpk_str = line.strip().split('\t')
		rpk = float(rpk_str)
		gene_fam_parts = gene_fam.split('|')
		if len(gene_fam_parts) == 1:
			pathway_data[gene_fam_parts[0]] = rpk
		elif len(gene_fam_parts) == 2:
			if gene_fam_parts[1] not in pathway_species_data:
				pathway_species_data[gene_fam_parts[1]] = {}
			pathway_species_data[gene_fam_parts[1]][gene_fam_parts[0]] = rpk

pathway_data_ordered = sorted(pathway_data.items(), key=lambda item: item[1])

species = 'g__Bacteroides.s__Bacteroides_uniformis_CAG_3'

specific_data_ordered = sorted(pathway_species_data[species].items(), key=lambda item: item[1], reverse=True)

for pathway, rpk in specific_data_ordered:
	print('%s: %s' % (rpk, pathway))
'''

fpath = '/u/home/d/daisyche/dbd/humann_data/hmp1-II_humann2_pathabundance-nrm-mtd-qcd.pcl.bz2'

hmp_samples = su.get_sample_names()

sample_pathway_data = {sample: {} for sample in hmp_samples}
sample_pathway_species_data = {sample: defaultdict() for sample in hmp_samples}

with bz2.BZ2File(fpath, 'r') as file:
	
	# Header lines
	header = file.readline()
	samples = header.strip().split('\t')[1:]
	
	hmp_sample_idxs = [samples.index(sample) for sample in hmp_samples]
	
	randsid_line = file.readline()
	visno_line = file.readline()
	areas = file.readline().strip().split('\t')[1:] # Gut, Oral, etc.
	sites = file.readline().strip().split('\t')[1:] # Stool, Tongue_dorsum, etc.
	snprnt_line = file.readline()
	gender_line = file.readline()
	wmsphase_line = file.readline()
	srs_line = file.readline()
	
	line_count = 0
	
	for line in file:
		if line_count % 100 == 0:
			print("On line %i..." % line_count)
		elems = line.strip().split('\t')
		pathway = elems[0]
		rpks = [float(elem) for elem in elems[1:]]
		pathway_parts = pathway.split('|')
		for i in hmp_sample_idxs:
			sample = samples[i]
			rpk = rpks[i]
			if len(pathway_parts) == 1:
				sample_pathway_data[sample][pathway_parts[0]] = rpk
			elif len(pathway_parts) == 2:
				if pathway_parts[1] not in sample_pathway_species_data[sample]:
					sample_pathway_species_data[sample][pathway_parts[1]] = {}
				sample_pathway_species_data[sample][pathway_parts[1]][pathway_parts[0]] = rpk
		line_count += 1

'''
pathway_data_ordered = sorted(sample_pathway_data.items(), key=lambda item: item[1])

species = 'g__Bacteroides.s__Bacteroides_uniformis_CAG_3'

specific_data_ordered = sorted(sample_pathway_species_data[sample][species].items(), key=lambda item: item[1], reverse=True)

for pathway, rpk in specific_data_ordered:
	print('%s: %s' % (rpk, pathway))
'''
	
# Pathway richness per sample per focal species
sample_focal_species_pathway_richness = {sample: defaultdict(int) for sample in hmp_samples}

# Pathway richness (set) per sample per focal species
sample_focal_species_pathway_richness_set = {sample: defaultdict(int) for sample in hmp_samples}

# Pathway richness of a focal species in a sample
sample_species_pathway_richness = {sample: defaultdict(int) for sample in hmp_samples}

# Pathway richness per sample (include duplicates across species)
sample_pathway_richness = defaultdict(int)

# Pathway richness per sample (set over all species)
sample_pathway_set_richness = {}

for sample in sample_pathway_species_data:
	
	pathways_set = set()
	focal_species_pathways_dict = defaultdict(list)
	
	for species in sample_pathway_species_data[sample]:
		focal_species_pathways = []
		for pathway in sample_pathway_species_data[sample][species]:
			rpk = sample_pathway_species_data[sample][species][pathway]
			if rpk > 0:
				sample_species_pathway_richness[sample][species] += 1
				sample_pathway_richness[sample] += 1
				pathways_set.add(pathway)
				focal_species_pathways_dict[species].append(pathway)
	
	for species in sample_pathway_species_data[sample]:
		sample_focal_species_pathway_richness[sample][species] = sample_pathway_richness[sample] - len(focal_species_pathways_dict[species])
		sample_focal_species_pathway_richness_set[sample][species] = len(pathways_set - set(focal_species_pathways_dict[species]))
	
	sample_pathway_set_richness[sample] = len(pathways_set)

# Store data
odir = '/u/home/d/daisyche/dbd/humann_data/'
output_fname = '%s/sample_pathway_richness.csv' % (odir)

output_file = open(output_fname, 'w')
output_file.write(','.join(['sample_id','pathway_richness','pathway_richness_set']) + '\n')
for sample in hmp_samples:
	items = [sample, sample_pathway_richness[sample], sample_pathway_set_richness[sample]]
	output_file.write(','.join(map(lambda x: str(x), items)) + '\n')

output_fname = '%s/sample_species_pathway_richness.csv' % (odir)

output_file = open(output_fname, 'w')
output_file.write(','.join(['sample_id','species','pathway_richness','sample_pathway_richness','sample_pathway_richness_set']) + '\n')
for sample in hmp_samples:
	for species in sample_species_pathway_richness[sample]:
		items = [sample, species, sample_species_pathway_richness[sample][species], sample_focal_species_pathway_richness[sample][species], sample_focal_species_pathway_richness_set[sample][species]]
		output_file.write(','.join(map(lambda x: str(x), items)) + '\n')

# Compare with alpha diversity
'''
pickle_dir = "%s/pickles" % config.data_directory
alpha_div_dict = pickle.load(open("%s/alpha_div_dict.pkl" % (pickle_dir), 'rb'))

samples = sample_pathway_richness.keys()
alpha_divs = [alpha_div_dict[sample] for sample in samples]
pathway_richnesses = [sample_pathway_richness[sample] for sample in samples]
pathway_richnesses_set = [sample_pathway_set_richness[sample] for sample in samples]

plt.plot(alpha_divs, pathway_richnesses, '.')
plt.xlabel("Sample Shannon alpha diversity")
plt.ylabel("Sample pathway richness (include duplicates across species)")
plt.show()

plt.close()

plt.plot(alpha_divs, pathway_richnesses_set, '.')
plt.xlabel("Sample Shannon alpha diversity")
plt.ylabel("Sample pathway richness (set)")
plt.show()
'''
