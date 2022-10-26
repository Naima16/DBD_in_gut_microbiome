genes_dir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/data/genes
species_list=$(cat $genes_dir/species_genes.txt)

gene_cov_file=$genes_dir/gene_coverage_per_sample_species.txt
gene_pres_file=$genes_dir/genes_present_per_sample_species.txt

touch $gene_cov_file
touch $gene_pres_file

# genes_summary.txt includes pangenome size and covered genes (at least 1 mapped read) info

IFS=$'\n'

for species in $species_list
do
	lines=$(cat $genes_dir/$species/genes_summary.txt | cut -f1-3 | tail -n +2)
	for line in $lines
	do
		echo -n $species >> $gene_cov_file
		echo -ne '\t' >> $gene_cov_file
		echo $line >> $gene_cov_file
	done
done

# genes_presabs.txt includes gene presence/absence info based on threshold on copy number, which is estimated by dividing average(?) read-depth of gene by median read-depth of 15 universal single copy genes

for species in $species_list
do
	ncols=$(bzcat $genes_dir/$species/genes_presabs.txt.bz2 | awk -F'\t' '{print NF; exit}')
	nrows=$(bzcat $genes_dir/$species/genes_presabs.txt.bz2 | wc -l)
	ngenes=$(($nrows - 1))
	
	for i in $(seq 2 $ncols) # Each column is a sample
	do		
		count=$(bzcat $genes_dir/$species/genes_presabs.txt.bz2 | cut -f $i | tail -n +2 | paste -sd+ | bc)
		
		echo -n $species >> $gene_pres_file
		echo -ne '\t' >> $gene_pres_file
		echo -n $ngenes >> $gene_pres_file
		echo -ne '\t' >> $gene_pres_file
		echo $count >> $gene_pres_file
	done
done

IFS=$' \t\n'
