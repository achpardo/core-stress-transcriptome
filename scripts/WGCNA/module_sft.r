##Picks soft thresholding power for module detection and coexpression network construction with WGCNA##
##Input files should be filtered with filter_data_wgcna.r first
##Input data should have genes in columns and samples in rows
##If data is input in .RData file then dataframe should be named "data_2"
##Follow up with module_detection.r
##--outputdir should have a "/" at the end

#Loads necessary libraries
library(optparse)
library(WGCNA)

#Important method for WGCNA to work
options(stringsAsFactors = FALSE)

#Initializes argments for command line
option_list = list(make_option(c("-i", "--input"), type="character", help="Input file"),
	make_option(c("-r", "--rdata"), help="Set to 'y' if input file is .RData"),
	make_option(c("-o", "--output"), help = "Name of output file - must be PDF"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
file = opt$input
out = opt$output
rdata = opt$rdata

if (rdata == 'n')
{
	#Reads in TPM data
	data <- read.table(file, sep = "\t", header = TRUE, row.names = 1)
	head(data)

	#Transposes data so WGCNA can read it
	data_2 = as.data.frame(t(data))
} else {

#Loads data
data = load(file)

}

#Determines proper soft thresholding power and plots results
powers = c(c(1:10), seq(from = 12, to=20, by=2))
print(powers)
sft = pickSoftThreshold(data_2, powerVector = powers, networkType = "signed", moreNetworkConcepts = TRUE, verbose = 5)
sizeGrWindow(9,5)
pdf(file = out, width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");

abline(h=0.80, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

print(warnings())