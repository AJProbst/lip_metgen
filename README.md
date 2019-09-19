# lip_metgen
Correlation and network analysis of metalipidomics and metagenomics data

code used for generating network analysis in https://www.biorxiv.org/content/10.1101/465690v1
usage:
organize two input files:
1. relative abundances of lipids, columns samples, rows lipid species
2. relative abundance of microbial species (e.g., from metagenomics), columns samples (same order as in file 1), rows lipid species

Workflow (use bash environment):
1. execute R script 01_correlation.R
2. execute ruby script 02_extract_clusters.rb
3. Visualization of output files is done in Cytoscape
