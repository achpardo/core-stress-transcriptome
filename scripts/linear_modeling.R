# using R version 4.2.1
# Purpose: Run a linear model to find out if the respective effects of genotype, treatment, and BioProject on PC1 of gene
## expression are significant (also include interaction effects, of course).
# Author: Anna Pardo
# Date initiated: Feb. 15, 2024

# load modules
library(RcppEigen)

# load data
expmd <- read.csv("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/principal_components_metadata.csv")

# try the linear model
#modelout = lm(PC1 ~ BioProject + Genotype + Treatment + Tissue + BioProject*Genotype + BioProject*Treatment + 
#                BioProject*Tissue + Genotype*Treatment + Genotype*Tissue + Treatment*Tissue + BioProject*Genotype*Treatment +
#                BioProject*Genotype*Tissue + BioProject*Treatment*Tissue + Genotype*Treatment*Tissue + 
#                BioProject*Genotype*Treatment*Tissue, data = expmd)

# try another model without the four-way interaction
## try one with just two-way interactions & first-order effects and see what you get

modelout <- lm(PC1 ~ Treatment + Tissue, data=expmd)

summary(modelout)

Anova(modelout)

# single-factor for genotype
gtmod <- lm(PC1 ~ Genotype, data=expmd)

summary(gtmod)
Anova(gtmod)

# try all four factors (no interactions)
fourfac <- lm(PC1 ~ Genotype + BioProject + Treatment + Tissue,data=expmd)
Anova(fourfac)
ffcoef <- as.data.frame(fourfac$coefficients)

# load dfgt20: same as expmd, but only for genotypes with >20 samples
dfgt20 <- read.csv("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/principal_components_metadata_gt20samp.csv")

# re-run four-factor no interaction model on this dataset
ffni <- lm(PC1 ~ Genotype + BioProject + Treatment + Tissue,data=dfgt20)
Anova(ffni)

# remove control for testing purposes
library(dplyr)
dfnocont <- dfgt20 %>% filter(Treatment != "Control")

ffni <- lm(PC1 ~ Genotype + BioProject + Treatment + Tissue,data=dfnocont)
Anova(ffni)
