import numpy as np
import bz2
import gzip

# DBD

fname='/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/data/snps/Bacteroides_vulgatus_57955/snps_depth.txt.bz2'

file = bz2.BZ2File(fname, 'r')

header1 = file.readline().split()

depths1_diff = []
depths1_same = []

for line in file:
	items = line.strip().split()
	depths1_diff.append(items[6])
	depths1_same.append(items[13])

a1 = np.array(depths1_diff[:330000]).astype('int')
b1 = np.array(depths1_same[:330000]).astype('int')

# DBD - intermediate file version

fname1 = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/debug/700013715/snps/output/Bacteroides_vulgatus_57955.snps.gz'

file = gzip.open(fname1, 'r')

header = file.readline()

depths1_diff_orig = [] # expect to be same as depths1_diff

for line in file:
	items = line.strip().split()
	depths1_diff_orig.append(items[5])

fname2 = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/debug/700014954/snps/output/Bacteroides_vulgatus_57955.snps.gz'

file = gzip.open(fname2, 'r')

header = file.readline()

depths1_same_orig = [] # expect to be same as depths1_same

for line in file:
	items = line.strip().split()
	depths1_same_orig.append(items[5])

# check if merge was correct
(np.array(depths1_same_orig[:1372450]) == np.array(depths1_same[:1372450])).sum() # should equal 1372450
(np.array(depths1_diff_orig[:1372450]) == np.array(depths1_diff[:1372450])).sum() # should equal 1372450

# yes, it was correct

# Mother-infant

fname='/u/project/ngarud/daisyche/old_mother_infant/data/snps/Bacteroides_vulgatus_57955/snps_depth.txt.bz2'

file = bz2.BZ2File(fname, 'r')

header2 = file.readline().split()

depths2_diff = []
depths2_same = []

for line in file:
	items = line.strip().split()
	depths2_diff.append(items[1])
	depths2_same.append(items[5])

a2 = np.array(depths2_diff[:330000]).astype('int')
b2 = np.array(depths2_same[:330000]).astype('int')

diffs_a = np.absolute(a2-a1)
diffs_b = np.absolute(b2-b1)

from matplotlib import pyplot as plt

plt.hist((a1-a2), bins=20)
plt.show()

plt.hist((b1-b2), bins=20)
plt.show()