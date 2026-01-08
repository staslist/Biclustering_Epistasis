# Biclustering_Epistasis

This is a companion repository for the EXTENSIVE GENETIC INTERACTIONS (EPISTASIS) LINKED TO ALCOHOL USE DISORDER IN A HIGH-RISK POPULATION manuscript.

Code related to (bi)clustering of SNP-SNP interaction results into interval-interval; gene-gene; pathway-pathway results.

BiClustering.py primarily contains code used in the bi-clustering algorithm. It also contains some code used to parse and annotate the interacting intervals output by bi-clustering.
Bio_Filter.py primarily contains code used to pre-filter the SNPs before conducting epistasis procedure.
hg19_to_hg38_ranges.py primarily contains code used to liftover the hg19 interacting interval loci to hg38.

For further explanation or detail please contact me via email or by posting an issue to this repository.

Based on: https://doi.org/10.1371/journal.pgen.1000782 Genome-Wide Association Data Reveal a Global Map of Genetic Interactions among Protein Complexes.

This code has been modified and further optimized for broader bioinformatics community use in the following GitHub repository: 
https://github.com/njdjyxz/geno-bicluster
