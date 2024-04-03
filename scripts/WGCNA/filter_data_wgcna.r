##Filters unecessary genes from gene expression (TPMs) file
##Input gene expression file should be tab delimited, with samples as columns, and gene names as rows
##Gene names should be the in the first column
##--outputdir should have "/" at the end

# install WGCNA
install.packages("WGCNA",repos="https://repo.miserver.it.umich.edu/cran/")

#Loads necessary libraries
library(optparse)
library(WGCNA)

#Important method for WGCNA to work
options(stringsAsFactors = FALSE)

#Initializes argments for command line
option_list = list(make_option(c("-i", "--input"), type="character", help="Gene expression input file"),
	make_option(c("-o", "--outputdir"), help = "Output directory"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
file = opt$input
out = opt$outputdir

setwd(out)

#Reads in TPM data
data <- read.table(file, sep = ",", header = TRUE, row.names = 1)
#print(head(data))
library(readr)
data <- read.table("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/BPcombat_TPM_coregenes_psyntissues_WGCNA.tsv")
rownames(data) <- data$Sample

#Transposes data so WGCNA can read it
data_2 = as.data.frame(t(data))
data_2 = t(data)
data_3 = data.frame(data_2)
print(dim(data_2))

#Filters data and writes new csv with only selected data

gsg = goodSamplesGenes(data, verbose = 1) #verbosity only refers to the actual output of the function - it does not effect how the function works!
gsg$allOK
if(!gsg$allOK)
{
	data_2 = data_2[gsg$goodSamples, gsg$goodGenes]
}
write.table(t(data_2), file = "data_used_for_network.txt", quote = FALSE, row.names = TRUE, sep = "\t")
print(dim(data_2))

#Hierarchal clustering to detect outliers
sampleTree = hclust(dist(data_2), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#Saves input data
save(data_2, file = "data_input.RData")
