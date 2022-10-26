Code for investigating the "diversity begets diversity" hypothesis in gut microbiota of healthy human adults, adapted from https://github.com/benjaminhgood/microbiome_evolution

Specialized for longitudinal data from Poyet

Overview of data:
- Source: https://www.ebi.ac.uk/ena/browser/view/PRJNA544527
- 5363 samples for whole PRJNA544527 study
- 7758 sample IDs in supplementary metadata table from https://www.nature.com/articles/s41591-019-0559-3
- 402 samples in Poyet_run_accessions.txt
- 401 samples we have MIDAS output for
- Only look at individuals ae, am, an, ao
- Individual ae (M, age 21): 59 samples
- Individual am (M, age 28): 206 samples (206 in Poyet_run_accessions.txt)
- Individual an (F, age 22): 63 samples
- Individual ao (M, age 37): 73 samples (74 in Poyet_run_accessions.txt)
	- The one missing sample is SRR9224118 in ao

To set up for use:
- Download a copy of the MIDAS database of reference genomes
- Download PATRIC kegg/feature files for all genomes in the MIDAS database
- Obtain MIDAS output for HMP1-2 dataset (should be a single directory containing merged genes, snps, species output)
- Run postprocessing scripts (postprocess/calculate_...) if necessary
- Make a symbolic link to utils under postprocess
