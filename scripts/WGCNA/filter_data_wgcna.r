##Filters unecessary genes from gene expression (TPMs) file
##Input gene expression file should be tab delimited, with samples as columns, and gene names as rows
##Gene names should be the in the first column
##--outputdir should have "/" at the end

# install WGCNA
#install.packages("WGCNA",repos="https://repo.miserver.it.umich.edu/cran/")

# install GO.db
#if (!require("BiocManager",quietly=TRUE))
#    install.packages("BiocManager")

#BiocManager::install("GO.db")
#library(GO.db)

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
data <- read.table(file, sep = "\t", header = TRUE, row.names = 1)

# read in TPM data - ondemand RStudio session
data <- read.table("/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/WGCNA_01-Apr-2024/BPcombat_TPM_allsamp_WGCNA.tsv",
                   sep = "\t", header = TRUE, row.names = 1)
#rownames(data) <- data$GeneID
#print(head(data))

#Transposes data so WGCNA can read it
data_2 = as.data.frame(t(data))
print(dim(data_2))

#Filters data and writes new csv with only selected data
gsg = goodSamplesGenes(data_2, verbose = 1) #verbosity only refers to the actual output of the function - it does not effect how the function works!
gsg$allOK
if(!gsg$allOK)
{
	data_2 = data_2[gsg$goodSamples, gsg$goodGenes]
}
write.table(t(data_2), file = "/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/WGCNA_01-Apr-2024/data_used_for_network.txt", quote = FALSE, row.names = TRUE, sep = "\t")
print(dim(data_2))
print(head(data_2))

#Hierarchal clustering to detect outliers
sampleTree = hclust(dist(data_2), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#Saves input data
save(data_2, file = "/mnt/research/VanBuren_Lab/02_users/Anna_Haber/core-stress-transcriptome/01_data/WGCNA_01-Apr-2024/data_input.RData")
