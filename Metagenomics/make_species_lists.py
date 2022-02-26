import bz2

file = bz2.BZ2File('../data/species/relative_abundance.txt.bz2', 'r')

header = file.readline()
samples = header.strip().split('\t')[1:]
sample_species_dict = {sample: [] for sample in samples}

for line in file:
	items = line.strip().split('\t')
	species = items[0]
	abs = [float(item) for item in items[1:]]
	for i in range(len(samples)):
		sample = samples[i]
		ab = abs[i]
		if ab > 0:
			sample_species_dict[sample].append(species)

for sample in sample_species_dict:
	output_file = open('../data/species_lists/%s_species_list.txt' % sample, 'w')
	species_list = sample_species_dict[sample]
	for species in species_list:
		output_file.write("%s\n" % species)
	print("Sample %s has %i species" % (sample, len(species_list)))
