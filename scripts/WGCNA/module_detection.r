##Performs hierarchal clustering on gene expression data using WGCNA
##Determine soft thresholding power using module_sft.r
##Input data loaded from .RData generated in module_sft.r

#Loads necessary libraries
library(optparse)

library(WGCNA)

#Important for WGNCA to work
options(stringsAsFactors = FALSE)

wd = getwd()

#Initializes argumnets for command line
option_list = list(make_option(c("-i", "--input"), help="Gene expression input as a .RData file"), 
	make_option(c("-s", "--sftpower"), type = "integer", help = "Number to indicate soft thresholding power"),
	make_option(c("-o", "--output"), help = "Output directory to store module data"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
file = opt$input
sft_power = opt$sftpower
out = opt$output

setwd(out)

#Reads input data obtained from module_sft.r
data = load("//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/WGCNA_output/data_input.RData")

#Defines soft thresholding power and creates adjacency matrix
softPower = 9
adjacency = adjacency(data_2, power = softPower, type = 'signed')

save(adjacency, file="adjacency_matrix.RData")

#Creates TOM matrix from adjacency
TOM = TOMsimilarity(adjacency)
dissTOM = TOMdist(adjacency)
print("hello")
#Performs hierarchal clustering and plots results
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
#plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

minModuleSize = 100
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#Clusters genes into modules
MEList = moduleEigengenes(data_2, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7,6)
pdf(file = "me_dendro.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
dev.off()

#Merges similar modules
merge = mergeCloseModules(data_2, dynamicColors, cutHeight = MEDissThres, verbose = 3)
head(merge)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

#Plots final dendrogram
sizeGrWindow(12, 9)
pdf(file = "gene_dendro.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

#Preps and saves module information
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
#save(MEs, moduleLabels, moduleColors, geneTree, TOM, dissTOM, file = "module_info.RData")
save(MEs, moduleLabels, moduleColors, geneTree, file = "module_info_no_TOM.RData")

setwd(wd)
