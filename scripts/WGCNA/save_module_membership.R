# Purpose: Get & save the module membership information.
# Author: Anna Pardo
# Date initiated: Apr. 6, 2024

# import modules
library(WGCNA)

# load module information
load("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_Apr2024/modules/module_info_no_TOM.RData")
# load datExpr
load("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_Apr2024/data_input.RData")

geneModuleMembership = as.data.frame(cor(data_2, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 1981))

datKME = signedKME(data_2, MEs, outputColumnName = "kME")
write.table(datKME,file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_Apr2024/module_membership_values.tsv",
            sep = "\t",row.names = TRUE, col.names = TRUE)
