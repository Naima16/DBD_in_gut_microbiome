# Investigating 'amount of evolution' vs alpha diversity
# For HMP adults only

import matplotlib	 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils import sample_utils as su, config
from utils.plot_utils import list_to_colors, colors_to_legend_elements

import sys, numpy as np
from collections import defaultdict

# ==========================================================
# Standard header to read in argument information
# ==========================================================
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--sweep-type', type=str, help="Full or partial sweep", default="full")

args = parser.parse_args()
sweep_type = args.sweep_type

if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep-type. Choose from full, partial")
# ==========================================================

# Load pickled data

import pickle
pickle_dir = "%s/pickles" % config.data_directory

alpha_div_dict = pickle.load(open("%s/alpha_div_dict.pkl" % (pickle_dir), 'rb'))
host_snp_change_dict = pickle.load(open("%s/host_snp_change_dict_%s.pkl" % (pickle_dir, sweep_type), 'rb'))
host_species_snp_change_dict = pickle.load(open("%s/host_species_snp_change_dict_%s.pkl" % (pickle_dir, sweep_type), 'rb'))
host_alpha_diversity_dict = pickle.load(open("%s/host_alpha_diversity_dict.pkl" % (pickle_dir), 'rb'))
host_change_type_dict = pickle.load(open("%s/host_change_type_dict_%s.pkl" % (pickle_dir, sweep_type), 'rb'))
sample_species_polymorphism_dict = pickle.load(open("%s/sample_species_polymorphism_dict_%s.pkl" % (pickle_dir, sweep_type), 'rb'))

# Set up proportion replacement vs modification for plot
# Threshold on at least 10 species present
prop_replace = []
prop_mod = []
prop_repmod_alpha_diversities = []

for subject in host_alpha_diversity_dict:
	
	if subject not in host_change_type_dict['none']:
		none_count = 0
	else:
		none_count = host_change_type_dict['none'][subject]
	
	if subject not in host_change_type_dict['replace']:
		replace_count = 0
	else:
		replace_count = host_change_type_dict['replace'][subject]
	
	if subject not in host_change_type_dict['mod']:
		mod_count = 0
	else:
		mod_count = host_change_type_dict['mod'][subject]
	
	total_count = replace_count + mod_count + none_count
	
	if total_count > 5:
		prop_replace.append(float(replace_count)/total_count)
		prop_mod.append(float(mod_count)/total_count)
		prop_repmod_alpha_diversities.append(host_alpha_diversity_dict[subject][0])

# Set up modification SNP change data for plot
host_alpha_diversities = []
host_snp_change_counts = []

for subject in host_alpha_diversity_dict:
	host_alpha_diversities.append(host_alpha_diversity_dict[subject][0])
	host_snp_change_counts.append(host_snp_change_dict[subject])

host_species_alpha_diversities = []
host_species_snp_change_counts = []
host_species_species = []

for subject, species in host_species_snp_change_dict:
	if subject in host_alpha_diversity_dict:
		host_species_alpha_diversities.append(host_alpha_diversity_dict[subject][0])
		host_species_snp_change_counts.append(host_species_snp_change_dict[(subject, species)])
		host_species_species.append(species)

# Set up within-sample polymorphism data for plot
sample_species_alpha_diversities = []
sample_species_polymorphisms = []

sample_alpha_diversities = []
sample_polymorphisms = []

for sample in sample_species_polymorphism_dict:
	if sample in alpha_div_dict:
		sample_alpha_diversities.append(alpha_div_dict[sample])
		agg_polymorphism = 0
		for species in sample_species_polymorphism_dict[sample]:
			sample_species_alpha_diversities.append(alpha_div_dict[sample])
			sample_species_polymorphisms.append(sample_species_polymorphism_dict[sample][species])
			agg_polymorphism += sample_species_polymorphism_dict[sample][species]
		sample_polymorphisms.append(agg_polymorphism/len(sample_species_polymorphism_dict[sample]))

# Plot 0: sample-species alpha diversity vs. polymorphism
fig, ax = plt.subplots()
ax.scatter(sample_species_alpha_diversities, sample_species_polymorphisms)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Within-sample polymorphism per sample-species pair")
ax.set_title("Alpha diversity vs polymorphism, HMP adults (n=%i)" % len(sample_species_polymorphisms))

fig.savefig('%s/alpha_div_vs_polymorphism_hmp.png' % (config.analysis_directory), bbox_inches='tight')

# Plot 0.5: sample alpha diversity vs. polymorphism (averaged across species)
fig, ax = plt.subplots()
ax.scatter(sample_alpha_diversities, sample_polymorphisms)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Within-sample polymorphism per sample, averaged across species")
ax.set_title("Alpha diversity vs polymorphism, HMP adults (n=%i)" % len(sample_polymorphisms))

fig.savefig('%s/alpha_div_vs_polymorphism_avg_hmp.png' % (config.analysis_directory), bbox_inches='tight')

# Plot 1: one point per host, aggregate across species
fig, ax = plt.subplots()
ax.scatter(host_alpha_diversities, host_snp_change_counts)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of SNP changes, aggregated across species")
ax.set_title("Alpha diversity vs number of SNP changes, HMP adults (n=%i)" % len(host_snp_change_counts))

fig.savefig('%s/alpha_div_vs_snp_changes_hmp_%s_agg.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 2: one point per host-species pair
fig, ax = plt.subplots()
colors = list_to_colors(host_species_species)
ax.scatter(host_species_alpha_diversities, host_species_snp_change_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of SNP changes per host-species pair")
ax.set_title("Alpha div. vs # SNP changes, HMP adults (n=%i)" % len(host_species_snp_change_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_snp_changes_hmp_%s_ind.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 3: Alpha diversity vs. fraction modification
fig, ax = plt.subplots()
ax.scatter(prop_repmod_alpha_diversities, prop_mod, color='r', label="Modification")
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Fraction of present species which are modified")
ax.set_title("Alpha diversity vs. fraction of modification events, HMP adults (n=%i)" % len(prop_mod))

fig.savefig('%s/alpha_div_vs_prop_mod_hmp_%s.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 4: Alpha diversity vs. fraction replacement
fig, ax = plt.subplots()
ax.scatter(prop_repmod_alpha_diversities, prop_replace, color='b', label="Replacement")
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Fraction of present species which are replaced")
ax.set_title("Alpha diversity vs. fraction of replacement events, HMP adults (n=%i)" % len(prop_replace))

fig.savefig('%s/alpha_div_vs_prop_replace_hmp_%s.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')
