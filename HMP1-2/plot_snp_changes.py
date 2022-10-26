from utils import config
import pickle, numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cmx
from math import exp, log

# Plotting utils
def calculate_unnormalized_survival_from_vector(counts):
	counts = sorted(counts)
	xs = []
	ns = []
	ns_cur = len(counts)
	min = -1
	for count in counts:
		if count > min:
			ns.append(ns_cur) # Number of elements greater or equal
			xs.append(count)
			min = count
		ns_cur -= 1
	xs.append(xs[len(xs)-1]+1)
	ns.append(0)
	return xs, numpy.array(ns)

# Load the pickled SNP change data

pooled_snp_pickle_fn = '%s/pickles/pooled_snp_change.pkl' % (config.data_directory)
pooled_between_snp_pickle_fn = '%s/pickles/pooled_between_snp_change.pkl' % (config.data_directory)

pooled_snp_change_distribution = pickle.load(open(pooled_snp_pickle_fn, 'rb'))
pooled_between_snp_change_distribution = pickle.load(open(pooled_between_snp_pickle_fn, 'rb'))

# Plot! (unrelated adult and adult-adult survival curves)
fig_snp, ax_snp = plt.subplots(figsize=(10,6))

colormap = cmx.get_cmap('jet', 2)
colors = [colormap(x) for x in numpy.array([x for x in range(0,2)])/2.0]

ax_snp.set_xscale('log')
ax_snp.set_yscale('log')
ax_snp.set_ylabel('Fraction comparisons $\geq n$', fontsize=11)
ax_snp.set_xlabel('# SNP changes', fontsize=11)

ax_snp.spines['top'].set_visible(False)
ax_snp.spines['right'].set_visible(False)
ax_snp.get_xaxis().tick_bottom()
ax_snp.get_yaxis().tick_left()

color_i = 0

# Within-host, adult
counts = []
for tp_pair in pooled_snp_change_distribution:
	counts += pooled_snp_change_distribution[tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = 'HMP: Adult 6 months' + (' (n=%d)' % ns[0])
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
ymin = 1.0/ns[0]
ymax = 1.3
ax_snp.set_ylim([ymin,ymax])
color_i += 1

# Unrelated adults
counts = []
for tp_pair in pooled_between_snp_change_distribution:
	counts += pooled_between_snp_change_distribution[tp_pair]

xs, ns = calculate_unnormalized_survival_from_vector(counts)
ax_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label="Unrelated adults", where='pre',zorder=4)
color_i += 1

ax_snp.legend(loc='best', frameon=True, fontsize=10, numpoints=1, ncol=1, handlelength=1)

modification_difference_threshold = config.modification_difference_threshold
replacement_difference_threshold = config.replacement_difference_threshold

# Now fill in the graphics
ax_snp.fill_between([1e-01,1], [ymin,ymin],[ymax,ymax],color='0.8',zorder=1)
ax_snp.fill_between([1e0,modification_difference_threshold],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
ax_snp.fill_between([replacement_difference_threshold,1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)

ax_snp.text( exp((log(1e05)+log(replacement_difference_threshold))/2), ymax*1.2, 'putative\nreplacement',fontsize=12,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
ax_snp.text( exp((log(1)+log(modification_difference_threshold))/2), ymax*1.2, 'putative\nmodification',fontsize=12,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)

fig_snp.savefig('%s/temporal_snp_changes_pooled.pdf' % (config.analysis_directory),bbox_inches='tight')