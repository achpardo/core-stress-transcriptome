# Purpose: Save module eigengene data from WGCNA.
# Author: Anna Pardo
# Date initiated: Apr. 10, 2024

library(WGCNA)

# load module information
load("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_Apr2024/modules/module_info_no_TOM.RData")

# save MEs dataframe
write.csv(MEs,"//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_Apr2024/module_eigengenes.csv",row.names = TRUE)
