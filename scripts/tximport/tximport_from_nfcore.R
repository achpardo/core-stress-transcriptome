# Title     : Tximport
# Objective : Map transcript level expression to gene level and calculate TPM (using output from nf-core pipeline)
# Created by: Jeremy (modified by Anna)
# Created on: 2021-04-27 (modified 2022-01-10 for sorghum multistress project; modified again 2022-02-14 for core stress project; 
## and again for use with nf-core rnaseq pipeline)
# deal with command line arguments
args = commandArgs(TRUE)
sampleTable = args[1]
# load required packages
library(readr)
library(dplyr)
library(tximport)
library(jsonlite)
output_dir <- "/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/03_processed_RNAseq_data"
# read in sample table containing salmon file paths
sampleTable <- read_csv(sampleTable)
head(sampleTable)
# extract name of BioProject for later use
bioproject <- unique(sampleTable$BioProject)
bioproject
#read in tx2gene file
tx2gene = read_delim("/mnt/scratch/haberan2/Core_Stress_Response/01_pipeline_outputs/May2_P39_and_others_nfcore/salmon/salmon_tx2gene.tsv",col_names= F, delim= "\t")
head(tx2gene)
# remove third column from tx2gene
tx2gene <- as.data.frame(tx2gene[,1:2])

# get vector of file paths
files = sampleTable$File
# run tximport
txi_ls = tximport(files = files,type="salmon", tx2gene=tx2gene,countsFromAbundance="lengthScaledTPM")
txi = tximport(files = files,type="salmon", tx2gene=tx2gene)
#get TPM
TPM_LS = as.data.frame(txi_ls$abundance)
colnames(TPM_LS) = sampleTable$Sample
TPM_LS$GeneID = rownames(TPM_LS)
Counts_LS = as.data.frame(txi_ls$counts)
colnames(Counts_LS) = sampleTable$Sample
Counts = as.data.frame(txi$counts)
colnames(Counts) = sampleTable$Sample
#save the results
save.image(paste0(output_dir,"/txi_",bioproject,"7-Jun-2023.RData"))
