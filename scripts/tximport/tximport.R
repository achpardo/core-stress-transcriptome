# Title     : Tximport
# Objective : Map transcript level expression to gene level and calculate TPM
# Created by: Jeremy (modified by Anna)
# Created on: 2021-04-27 (modified 2022-01-10 for sorghum multistress project; modified again 2022-02-14 for core stress project [rice])
# load required packages
library(readr)
library(dplyr)
library(tximport)
library(jsonlite)
output_dir <- "/mnt/gs21/scratch/haberan2/Core_Stress_Response/01_pipeline_outputs/B73_HS_Mar28/alltximport"
# Read in sample table containing salmon filepaths
sampleTable <- read_csv("/mnt/research/VanBuren_Lab/02_users/Anna_Haber/02_core_stress/01_maize/B73HS_salmonmap.csv")
head(sampleTable)
#read in tx2gene file
tx2gene = read_delim("/mnt/gs21/scratch/haberan2/Core_Stress_Response/01_pipeline_outputs/B73_HS_Mar28/batch3/tx2gene.txt",col_names= F, delim= "\t")
head(tx2gene)

# get vector of file paths
files = sampleTable$File
# run tximport
txi_ls = tximport(files = files,type="salmon", tx2gene=tx2gene,countsFromAbundance="lengthScaledTPM")
txi = tximport(files = files,type="salmon", tx2gene=tx2gene)
#get TPM
TPM_LS = as.data.frame(txi_ls$abundance)
colnames(TPM_LS) = sampleTable$Sample
TPM_LS$GeneID = rownames(TPM_LS)
write_delim(TPM_LS,paste0(output_dir,"/TPM_LS_B73HS_29-Mar-2023.txt"),delim="\t")
Counts_LS = as.data.frame(txi_ls$counts)
colnames(Counts_LS) = sampleTable$Sample
Counts = as.data.frame(txi$counts)
colnames(Counts) = sampleTable$Sample
#save the results
save.image(paste0(output_dir,"/txi_B73HS_29-Mar-2023.RData"))
