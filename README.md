## Community diversity is associated with intra-species genetic diversity and gene loss in the human gut microbiome (Madi et al. 2022)

Code for metagenomic pipeline (MIDAS) and post-Midas filters, statistics and figures in Madi et al. 2022. 
https://doi.org/10.1101/2022.03.08.483496

Two publicly available metagenomic datasets were examined in this study:
1. Human Microbiome Project Consortium 2012 and Lloyd-Price et al. (2017)
(URL: https://aws.amazon.com/datasets/human-microbiome-project/)
2. Poyet et al. 2019 (NCBI accession number PRJNA544527)

Both datasets were initially processed using the software MIDAS (Nayfach et al. 2016, https://github.com/snayfach/MIDAS) for estimating species, gene and SNV content of metagenomic samples.

The `Metagenomics` subdirectories contain code for postprocessing MIDAS output (adapted from Garud and Good et al. 2019, https://github.com/benjaminhgood/microbiome_evolution). The `Statistics` folder contains R scripts used for statistical analyses and figure plotting.

