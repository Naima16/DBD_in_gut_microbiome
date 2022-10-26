# Investigating gene gain/loss vs alpha diversity
# For HMP adults only

import matplotlib	 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils import sample_utils as su, config
from utils.plot_utils import list_to_colors, colors_to_legend_elements

import sys, numpy as np
from collections import defaultdict

# Load pickled data

import pickle
pickle_dir = "%s/pickles" % config.data_directory

alpha_div_dict = pickle.load(open("%s/alpha_div_dict.pkl" % (pickle_dir), 'rb'))
host_alpha_diversity_dict = pickle.load(open("%s/host_alpha_diversity_dict.pkl" % (pickle_dir), 'rb'))
host_gene_gain_dict = pickle.load(open("%s/host_gene_gain_dict.pkl" % (pickle_dir), 'rb'))
host_gene_loss_dict = pickle.load(open("%s/host_gene_loss_dict.pkl" % (pickle_dir), 'rb'))
host_species_gene_gain_dict = pickle.load(open("%s/host_species_gene_gain_dict.pkl" % (pickle_dir), 'rb'))
host_species_gene_loss_dict = pickle.load(open("%s/host_species_gene_loss_dict.pkl" % (pickle_dir), 'rb'))

# Set up gene gain/loss data for plot

host_alpha_diversities = []
host_gene_gain_counts = []
host_gene_loss_counts = []

host_alpha_div_diffs = []
host_gene_gain_counts_for_diffs = []
host_gene_loss_counts_for_diffs = []

for subject in host_alpha_diversity_dict:
	alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
	gain_count = host_gene_gain_dict[subject]
	loss_count = host_gene_loss_dict[subject]
	
	if alpha_div1 != None:
		host_alpha_diversities.append(alpha_div1)
		host_gene_gain_counts.append(gain_count)
		host_gene_loss_counts.append(loss_count)
	
	if alpha_div1 != None and alpha_div2 != None:
		alpha_div_diff = alpha_div2 - alpha_div1
		host_alpha_div_diffs.append(alpha_div_diff)
		host_gene_gain_counts_for_diffs.append(gain_count)
		host_gene_loss_counts_for_diffs.append(loss_count)

host_gene_change_counts = np.array(host_gene_gain_counts) + np.array(host_gene_loss_counts)
host_gene_change_counts_for_diffs = np.array(host_gene_gain_counts_for_diffs) + np.array(host_gene_loss_counts_for_diffs)

host_species_gene_gain_alpha_divs = []
host_species_gene_gain_counts = []
host_species_gene_gain_species = []

host_species_gene_gain_alpha_div_diffs = []
host_species_gene_gain_counts_for_diffs = []
host_species_gene_gain_species_for_diffs = []

host_species_gene_loss_alpha_divs = []
host_species_gene_loss_counts = []
host_species_gene_loss_species = []

host_species_gene_loss_alpha_div_diffs = []
host_species_gene_loss_counts_for_diffs = []
host_species_gene_loss_species_for_diffs = []

for subject, species in host_species_gene_gain_dict:
	
	gain_count = host_species_gene_gain_dict[(subject, species)]
	
	if subject in host_alpha_diversity_dict:
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
		
		if alpha_div1 != None:
			host_species_gene_gain_alpha_divs.append(alpha_div1)
			host_species_gene_gain_counts.append(gain_count)
			host_species_gene_gain_species.append(species)
		
		if alpha_div1 != None and alpha_div2 != None:
			alpha_div_diff = alpha_div2 - alpha_div1
			host_species_gene_gain_alpha_div_diffs.append(alpha_div_diff)
			host_species_gene_gain_counts_for_diffs.append(gain_count)
			host_species_gene_gain_species_for_diffs.append(species)

for subject, species in host_species_gene_loss_dict:
	
	loss_count = host_species_gene_loss_dict[(subject, species)]
	
	if subject in host_alpha_diversity_dict:
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
		
		if alpha_div1 != None:
			host_species_gene_loss_alpha_divs.append(alpha_div1)
			host_species_gene_loss_counts.append(loss_count)
			host_species_gene_loss_species.append(species)
		
		if alpha_div1 != None and alpha_div2 != None:
			alpha_div_diff = alpha_div2 - alpha_div1
			host_species_gene_loss_alpha_div_diffs.append(alpha_div_diff)
			host_species_gene_loss_counts_for_diffs.append(loss_count)
			host_species_gene_loss_species_for_diffs.append(species)

# Plot 1: one point per host, aggregate across species
fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_gain_counts, 'ro', label='Gene gains')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene gains, aggregated across species")
ax.set_title("Alpha diversity vs number of gene gains, HMP adults (%i subjects)" % len(host_gene_gain_counts))

fig.savefig('%s/alpha_div_vs_gene_gains_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_loss_counts, 'bo', label='Gene losses')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene losses, aggregated across species")
ax.set_title("Alpha diversity vs number of gene losses, HMP adults (%i subjects)" % len(host_gene_loss_counts))
fig.savefig('%s/alpha_div_vs_gene_losses_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_change_counts, 'go', label='Combined')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene changes (gain or loss), aggregated across species")
ax.set_title("Alpha diversity vs number of gene changes, HMP adults (%i subjects)" % len(host_gene_change_counts))

fig.savefig('%s/alpha_div_vs_gene_changes_hmp_agg.png' % (config.analysis_directory))

# =================================================================
# Plot 2: one point per host-species pair
# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_gain_species)
ax.scatter(host_species_gene_gain_alpha_divs, host_species_gene_gain_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene gains per host-species pair")
ax.set_title("Alpha div. vs # gene gains,\nHMP adults (n=%i)" % len(host_species_gene_gain_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_gain_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_gene_gains_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_loss_species)
ax.scatter(host_species_gene_loss_alpha_divs, host_species_gene_loss_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene losses per host-species pair")
ax.set_title("Alpha div. vs # gene losses,\nHMP adults (n=%i)" % len(host_species_gene_loss_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_loss_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_gene_loss_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================
# Plot 3: one point per host, alpha div diffs
# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_gain_counts_for_diffs, 'ro', label='Gene gains')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene gains, aggregated across species")
ax.set_title("Alpha diversity change vs # gene gains,\nHMP adults (%i subjects)" % len(host_gene_gain_counts_for_diffs))

plt.grid(True,which="both", linestyle='--')
fig.savefig('%s/alpha_div_diff_vs_gene_gains_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_loss_counts_for_diffs, 'bo', label='Gene losses')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene losses, aggregated across species")
ax.set_title("Alpha diversity change vs # gene losses,\nHMP adults (%i subjects)" % len(host_gene_loss_counts_for_diffs))

fig.savefig('%s/alpha_div_diff_vs_gene_losses_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_change_counts_for_diffs, 'go', label='Combined')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene changes (gain or loss), aggregated across species")
ax.set_title("Alpha diversity change vs # gene changes,\nHMP adults (%i subjects)" % len(host_gene_change_counts_for_diffs))

fig.savefig('%s/alpha_div_diff_vs_gene_changes_hmp_agg.png' % (config.analysis_directory))

# =================================================================
# Plot 4: one point per host-species pair, alpha div diffs
# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_gain_species_for_diffs)
ax.scatter(host_species_gene_gain_alpha_div_diffs, host_species_gene_gain_counts_for_diffs, c=colors, edgecolors='none')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene gains per host-species pair")
ax.set_title("Alpha div. vs # gene gains,\nHMP adults (%i host-species pairs)" % len(host_species_gene_gain_counts_for_diffs))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_gain_species_for_diffs), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_diff_vs_gene_gains_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_loss_species_for_diffs)
ax.scatter(host_species_gene_loss_alpha_div_diffs, host_species_gene_loss_counts_for_diffs, c=colors, edgecolors='none')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene losses per host-species pair")
ax.set_title("Alpha div. vs # gene losses,\nHMP adults (%i host-species pairs)" % len(host_species_gene_loss_counts_for_diffs))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_loss_species_for_diffs), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_diff_vs_gene_loss_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')