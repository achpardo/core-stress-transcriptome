# using R version 4.2.1
# Purpose: Run a linear model to find out if the respective effects of genotype, treatment, and BioProject on PC1 of gene
## expression are significant (also include interaction effects, of course).
# Author: Anna Pardo
# Date initiated: Feb. 15, 2024

# load modules
library(dplyr)

# load data
expmd <- read.csv("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/principal_components_metadata.csv")

# try the linear model
modelout = lm(PC1 ~ BioProject + Genotype + Treatment + Tissue + BioProject*Genotype + BioProject*Treatment + 
                BioProject*Tissue + Genotype*Treatment + Genotype*Tissue + Treatment*Tissue + BioProject*Genotype*Treatment +
                BioProject*Genotype*Tissue + BioProject*Treatment*Tissue + Genotype*Treatment*Tissue + 
                BioProject*Genotype*Treatment*Tissue, data = expmd)
