#!/usr/bin/env python3
from StrainFinder import *
import numpy as np
import pandas as pd
import pickle
import os
import sys
import time

def get_best_EM(sample,species):
    
    #fns = [f"{stem}/{sample}/{species}/{species}{n}"  for n in range(1,5)]
    fns = [f"{stem}/{sample}/{species}/{file}" for file in os.listdir(f"{stem}/{sample}/{species}/")]
    # Load EM objects
    ems = [pickle.load(open(fn, 'rb'), encoding='latin1') for fn in fns]

    # Get the best AIC in each EM object
    aics = [em.select_best_estimates(1)[0].aic for em in ems]

    # Select EM with the minimum AIC
    best_em = ems[np.argmin(aics)]
    
    return(best_em.select_best_estimates(1)[0])

stem="/u/scratch/r/rwolff/strainFinder/host_sf_freqs"
species=str(sys.argv[1])
samples = []
for sample in os.listdir(stem):
    if species in os.listdir(f"{stem}/{sample}"):
        samples.append(sample)

i = 1

strain_list = []
sample_list = []

for sample in samples:
    try:
        best = get_best_EM(sample,species)
        for i in range(1,len(best.p)+1):
            sample_list.append(sample)
            strain_list.append(i)
        print(sample)
    except:
        pass
    
outfile = pd.DataFrame(strain_list,index=sample_list,columns=["strain"])
outfile.to_csv(f"/u/home/r/rwolff/oligocolonization/StrainFinder/scripts/strain_number_files/{species}")