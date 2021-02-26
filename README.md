# Potamilus Ecology and Evolution 
Molecular data for Smith et al. 2021 - Establishing conservation units to promote recovery of two threatened freshwater mussel species (Bivalvia: Unionida: Potamilus)

Information about files:
**amponly_ind_assignments.csv** - Designations to create the dataframe for P. amphichaenus only 
**Pot_cov.csv** - covfile to link metadata to genotype data
**Potamilus_RAD.R** - R script used for data upload, filtering, and analyses. Datasets for external analyses were created using this script (e.g., fasta, structure)  
**Report_DPota20-4942_SNP_2.csv** - SNP calls from DArTSoft14 pipeline and associated metadata. This file contains 65,465 SNPs for 91 individuals (after exluding one poor sample) that were further filtered for analysis
**rm_PohiBra036.csv** - Cull sample PohiBra036 due to poor genotyping
**stronly_ind_assignments.csv** - Designations to create the dataframe for P. streckersoni only 

Notes:
- Raw reads are deposited in the SRA (BioProject ID: PRJNA663379) and sample accession numbers can be found in Table 1 of Smith et al. 2021
- ND1 sequences used for mitochondrial data analysis are deposited on GenBank and accession numbers can be found in Table 1 of Smith et al. 2021
- Data for mapping has to be downloaded from public repositories (USGS, Natural Earth)
