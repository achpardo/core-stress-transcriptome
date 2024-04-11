# Purpose: Get intramodular connectivity for each coexpression module & identify hub genes.
# Author: Anna Pardo
# Date initiated: Apr. 5, 2024

# load modules
library(WGCNA)

# load adjacency matrix
amat <- load("/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/WGCNA_01-Apr-2024/modules/adjacency_matrix.RData")

# load module information
mi <- load("/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/WGCNA_01-Apr-2024/modules/module_info_no_TOM.RData")

intconn <- intramodularConnectivity(adjMat = adjacency, colors = moduleColors)