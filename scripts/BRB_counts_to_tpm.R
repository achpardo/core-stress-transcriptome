# Purpose: Convert counts table from Brandon Webster (BRB-seq of low nitrogen stress) to TPM.
# Author: Anna Pardo
# Date initiated: May 11, 2023

# load modules
library(edgeR)
library(dplyr)
library(readr)
library(DGEobj.utils)
library(janitor)
library(purrr)

# load data
## first: load gene lengths matrix
geneLengths <- read_delim("../data/gene_lengths_B73V5.txt",delim = "\t")
head(geneLengths)

## second, load counts matrix
## this was acquired from Brandon Webster (Dr. Addie Thompson's lab, MSU) for his low nitrogen experiment
counts <- read_csv("../data/All_DGE.csv")
str(counts)

# data wrangling of counts matrix: need to generate unique sample names
# most of the sample names are unique, but one, "B73xPHZ51_LII_3", is duplicated.
# does this extend to the whole rows? subset counts df to only contain data for this sample.
counts %>% filter(Sample=="B73xPHZ51_LII_3")
# no, it is in two different BRBPools. 
# Solution: concatenate BRBPool to Sample.
counts$Sample <- paste0(counts$Sample,"_BRB",counts$BRBPool)

# slice off all metadata columns except sample
countsSub <- counts %>% select(-c(Plot,Range,Pass,Block,Pedigree,N_ResponseType,Treatments,DevStage,TreatBlocks,BRBPool))

# convert data from wide to long with the samples as the headers
countsSubT <- t(countsSub)
countsSubT <- row_to_names(dat=countsSubT,row_number = 1)

# extract row names as a column
countsSubT2 <- data.frame(countsSubT)
countsSubT2 <- tibble::rownames_to_column(countsSubT2,"GeneID")
head(countsSubT2)

# merge countsSubT2 with geneLengths
dfMerged <- merge(geneLengths,countsSubT2,by = "GeneID")
head(dfMerged)

# now split off the counts matrix for convertCounts() (all columns must be numeric)
countsMatrix <- dfMerged %>% select(-c(GeneID,Length))
str(countsMatrix)
# convert all columns to numeric
countsMatrixNum <- countsMatrix %>% mutate_all(function(x) as.numeric(as.character(x)))

# and split off a new geneLengths object
geneLengths2 <- dfMerged %>% select(c(GeneID,Length))

# all this data splitting just ensures that the dataframes are both in the same order.
# it's time to generate TPM using convertCounts()

brbtpm <- convertCounts(countsMatrix = as.matrix(countsMatrixNum),unit = "tpm",geneLength = geneLengths2$Length,log = FALSE,normalize = "none")
head(brbtpm)

# add back gene IDs and save as txt file
brbtpmDF <- data.frame(brbtpm)
brbtpmDF$GeneID <- dfMerged$GeneID
write_delim(brbtpmDF,file = "../data/Webster_TPM.txt",delim = "\t")
# also save the metadata from counts to a dataframe
meta <- counts %>% select(c(Plot,Range,Pass,Block,Pedigree,N_ResponseType,Treatments,DevStage,TreatBlocks,BRBPool,Sample))
write_delim(meta,file = "../data/Webster_metadata.txt",delim = "\t")
