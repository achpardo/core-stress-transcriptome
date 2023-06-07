# Purpose: Run DESeq2 (differential expression) on samples for core stress 
## project.
# Author: Anna Pardo
# Date initiated: June 5, 2023

# load required packages
library(readr)
library(tximport)
library(DESeq2)
library(dplyr)

# set working directory
setwd("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/scripts/differential_expression/")

###### Run DESeq2 from tximport output #######

# load tximport RData file
load("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/txi_forDE_06-Jun-2023.RData")

# data wrangling
## load in list of samples for DE
samp <- read_delim("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/samples_for_de.txt",col_names = F,delim = "\t")
samp <- samp$X1

## load in metadata
md <- read_csv("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/metadata_for_DESeq_samples.csv",col_names = T)
head(md)

# subset metadata to only the samples found in the sample list
mdsub <- md[which(md$Sample %in% samp),]

# generate unique sample condition identifiers containing the following information:
#   BioProject
#   Genotype
#   Treatment
#   Duration_hours (time point)
#   Concentration_mM
#   Concentration
#   Developmental_stage
#   Tissue
mdsub$identifier <- paste(mdsub$BioProject,mdsub$Genotype,mdsub$Treatment,mdsub$Duration_hours,mdsub$Concentration_mM,mdsub$Concentration,
                          mdsub$Developmental_stage,mdsub$Tissue,sep = "_")

# run DESeq
DESeq_obj <- DESeqDataSetFromTximport(txi = txi,colData = mdsub,design = ~identifier)
dds_obj <- DESeq(DESeq_obj)
