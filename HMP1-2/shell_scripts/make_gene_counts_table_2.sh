genes_dir=/u/home/d/daisyche/dbd/data/genes
gene_counts_dir=/u/home/d/daisyche/dbd/data/gene_counts/gene_counts_v2

species_genes=`cat $genes_dir/species_genes.txt`

output_file=/u/home/d/daisyche/dbd/data/gene_counts/gene_counts_per_sample_species_v2.txt

species=Fusobacterium_mortiferum_61584 # Random one
file=$gene_counts_dir/gene_counts_$species.csv
cat $file | head -1 >> $output_file

for species in $species_genes # Each column is a sample
do
	echo Working on species $species
	file=$gene_counts_dir/gene_counts_$species.csv
	cat $file | tail -n +2 >> $output_file
done
