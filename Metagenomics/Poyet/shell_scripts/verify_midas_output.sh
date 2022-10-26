#!/bin/bash

# ===================================================================================
# Verifies that all relevant samples have complete MIDAS intermediate output
# 
# Specifically, this script checks that:
# 	- there exist a genes, snps and species folder
# 	- under species, there exist log.txt, readme.txt, species_profile.txt
# 	- under genes, there exist log.txt, output, readme.txt, species.txt, summary.txt
# 	- there is no temp under genes
# 	- number of files in genes/output matches length of species.txt
# 	- under snps, there exist log.txt, output, readme.txt, species.txt, summary.txt
# 	- there is no temp under snps
# 	- number of files in species/output matches length of species.txt
# ===================================================================================

proj_dir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2

samples_file=$proj_dir/HMP1-2_samples_final.txt
samples=$(cat $samples_file)

midas_dir=$proj_dir/midas_output

snps_files='log.txt readme.txt species.txt summary.txt'
genes_files='log.txt readme.txt species.txt summary.txt'

for sample in $samples; do
	
	species_dir=$midas_dir/$sample/species
	snps_dir=$midas_dir/$sample/snps
	genes_dir=$midas_dir/$sample/genes
	
	# Check for species, snps, genes
	
	if [ ! -d $species_dir ] || [ ! -d $snps_dir ] || [ ! -d $genes_dir ];
	then echo $sample is missing either species, snps or genes; break; fi
	
	# Check species for log.txt, readme.txt, species_profile.txt
	
	if [ ! -f $species_dir/log.txt ] || [ ! -f $species_dir/readme.txt ] || [ ! -f $species_dir/species_profile.txt ];
	then echo $sample is missing species information;	fi
	
	# Check snps for log.txt, readme.txt, species.txt, summary.txt
	
	for file in $snps_files; do
		if [ ! -f $snps_dir/$file ];
		then echo $sample is missing snps information; break; fi
	done
	
	# Check snps for temp directory
	
	if [ -d $snps_dir/temp ];
	then echo $sample has temp directory in snps; fi
	
	# Check snps output for correct length and emptiness
	
	output_len=$(ls $snps_dir/output | wc -w)
	expect_len=$(cat $snps_dir/species.txt | wc -w)
	if [ ! $output_len -eq $expect_len ] || [ $output_len -eq 0 ];
	then echo $sample snps output is empty or the wrong length; fi
	
	# Check genes for log.txt, readme.txt, species.txt, summary.txt
	
	for file in $genes_files; do
		if [ ! -f $genes_dir/$file ];
		then echo $sample is missing genes information; break; fi
	done
	
	# Check genes for temp directory
	
	if [ -d $genes_dir/temp ];
	then echo $sample has temp directory in genes; fi
	
	# Check genes output for correct length and emptiness
	
	output_len=$(ls $genes_dir/output | wc -w)
	expect_len=$(cat $genes_dir/species.txt | wc -w)
	if [ ! $output_len -eq $expect_len ] || [ $output_len -eq 0 ];
	then echo $sample genes output is empty or the wrong length; fi
	
done
