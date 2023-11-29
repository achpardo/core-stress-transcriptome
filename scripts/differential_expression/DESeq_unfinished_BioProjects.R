# Purpose: Run DESeq2 (differential expression) on samples for core stress 
## project. Specifically, this script is for the following BioProjects that
## previously did not finish:
#   PRJNA414300
#   PRJNA420600
#   PRJNA611589
#   PRJNA687609
#   PRJNA877073
# Author: Anna Pardo
# Date initiated: June 5, 2023

# load required packages
library(readr)
library(tximport)
library(DESeq2)
library(dplyr)
library(stringr)

###### Run DESeq2 from tximport output #######

# load tximport RData file
load("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/txi_PRJNA6876097-Jun-2023.RData")

# data wrangling
## load in metadata
md <- read_csv("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/metadata_for_DESeq_samples.csv",col_names = T)
head(md)
md$identifier <- paste(md$Genotype,md$Treatment,md$Duration_hours,md$Concentration_mM,md$Concentration,md$Developmental_stage,md$Tissue,sep = "_")
# save md for use in python
write_csv(md,file = "//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/metadata_ID_DE.csv")

# subset metadata to only contain the BioProject of relevance
mdb <- md %>% filter(BioProject == "PRJNA687609")
head(mdb)
# for PRJNA414300 only: change "-" in "V1-V3" to "to"
mdb$Developmental_stage <- sub("-","to",mdb$Developmental_stage)
# for PRJNA611589 and PRJNA687609 only: remove spaces in genotype names and the dash in one of the genotypes
mdb$Genotype <- sub(" ","",mdb$Genotype)
mdb$Genotype <- sub("-","dash",mdb$Genotype)

# generate unique sample condition identifiers containing the following information:
#   Genotype
#   Treatment
#   Duration_hours (time point)
#   Concentration_mM
#   Concentration
#   Developmental_stage
#   Tissue
mdb$identifier <- paste(mdb$Genotype,mdb$Treatment,mdb$Duration_hours,mdb$Concentration_mM,mdb$Concentration,mdb$Developmental_stage,mdb$Tissue,sep = "_")

# run DESeq
DESeq_obj <- DESeqDataSetFromTximport(txi = txi,colData = mdb,design = ~identifier)
dds_obj <- DESeq(DESeq_obj)

get_sig_df = function(df,cont){
    df = mutate(df,GeneID = rownames(df),Contrast = cont)
    df = df[which(df$padj<0.05),]
    return(df)
}

expConditions <- unique(mdb$identifier)

# set up reference levels: anything that contains Control is a reference level
refLevels = list()
x=1

for(i in 1:length(expConditions)){
	if(str_detect(expConditions[i],"Control")==TRUE){
		refLevels[x] <- expConditions[i]
		x=x+1
	}
}

# the following code is mostly taken from Jeremy Pardo's DESeq2 script
# https://github.com/pardojer23/RNAseqV2/blob/master/RNAseqV2/R_Scripts/DESeq2.r
print(paste("Collecting results using the following treatment(s) as reference:",refLevels,sep=" "))

compare_list = list()
x=1
for (i in 1:length(refLevels)){
    for (j in 1:length(expConditions)){
        if (expConditions[j] != refLevels[i]){
            compare_list[x] = paste0(as.character(refLevels[i]),"-",as.character(expConditions[j]))
            x = x+1
        }
    }
}

df_list = list()
for (i in 1:length(compare_list)){
    print(compare_list[i][1])
    condition_list = unlist(strsplit(as.character(compare_list[i][1]),"-"))
    contrast = paste0(condition_list[1],"-V-",condition_list[2])
    print(contrast)
    df_list[[i]] = assign(paste0(contrast,"_sig"),get_sig_df(assign(contrast,as.data.frame(results(dds_obj,contrast=c("identifier",condition_list[2],condition_list[1]),alpha=0.05,pAdjustMethod = "fdr"))),contrast))
}
DE_Gene_df = bind_rows(df_list)

# set filename
filename = paste0("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/","PRJNA687609","_DEG_df.txt")
#write file
write_delim(DE_Gene_df,filename,delim="\t")
