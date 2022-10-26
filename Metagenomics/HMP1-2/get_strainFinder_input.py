import numpy as np
import pandas as pd
import sys
import os
import pickle

def get_counts(counts,depth_thresh):
    
    emp = []
    
    for i in range(counts.shape[0]):
        
        p = counts[i].split(',')
        elem = [int(e) for e in p]
        
        s = sum(e_s)

        ## here, we look for sites that have more than 0 reads mapping 
        ## to the second most common nucleotide e.g. polymorphic sites
        ## we append such sites to the list emp
        e_s = sorted(elem)
        if e_s[-2] > 0 and s > 20:
            emp.append(elem)

    y = np.array(emp)
    z = np.array([y])
    return z

cohort = "HMP1-2"
stem = f"/u/project/ngarud/Garud_lab/metagenomic_fastq_files/{cohort}/midas_output/"
num = int(sys.argv[1])-1

samples = os.listdir(stem)
sample = samples[num]

if not os.path.isdir(f"/u/scratch/r/rwolff/strainFinder/host_sf_inputs/{cohort}/{sample}"):
    os.makedirs(f"/u/scratch/r/rwolff/strainFinder/host_sf_inputs/{cohort}/{sample}",exist_ok=True)

# optionally, threshold on species prevalence
#species_file = pd.read_csv(f"/u/project/ngarud/Garud_lab/metagenomic_fastq_files/{cohort}/data/species/species_prevalence.txt.bz2",index_col=0,sep="\t",engine="python")
#species_file = species_file[species_file["prevalence"] >= 5]
species_file = [e.split(".")[0] for e in os.listdir(f"{stem}{sample}/snps/output/")]
species_list = list(species_file)
print(species_list)

for species in species_list:
    
    snps_loc = f"{stem}{sample}/snps/output/{species}.snps.gz"
    df = pd.read_csv(snps_loc,sep='\t')
    counts = np.array(df['count_atcg'])
    atcg = get_counts(counts,0)
    
    with open(f"/u/scratch/r/rwolff/strainFinder/host_sf_inputs/{cohort}/{sample}/{species}.pkl",'wb') as f:
        pickle.dump(atcg, f,protocol=2)
    
    sys.stdout.write("{0} finished! \n".format(species))
    sys.stdout.flush()
    
sys.stdout.write("{0} finished! \n \n \n".format(sample))
sys.stdout.flush()

