##### Exploratory section: figuring things out #######
library(limma)
z <- getGeneKEGGLinks("zma")
head(z)

library(OmnipathR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")

library(KEGGREST)
head(keggConv("zma","uniprot"))

library(biomaRt)

ensembl_plants = useMart("plants_mart",host="plants.ensembl.org")
datasets = listDatasets(ensembl_plants)
maize_data = useDataset("zmays_eg_gene",mart=ensembl_plants)
maize_mart = useMart("plants_mart",host="plants.ensembl.org", dataset = "zmays_eg_gene")

genes = getBM(attributes = c("ensembl_gene_id","entrezgene_id"),mart=maize_mart)
atr = listAttributes(maize_mart)
atr$name[37]

test_df = as.data.frame(genes)

########## Actual code: Pulling together NCBI & MaizeGDB gene IDs ########
ensembl_plants = useMart("plants_mart",host="plants.ensembl.org")
maize_data = useDataset("zmays_eg_gene",mart=ensembl_plants)
maize_mart = useMart("plants_mart",host="plants.ensembl.org", dataset = "zmays_eg_gene")
genes = getBM(attributes = c("ensembl_gene_id","entrezgene_id"),mart=maize_mart)
genes_df = as.data.frame(genes)

# drop NA rows in genes_df
library(dplyr)
library(tidyr)
gdf <- genes_df %>% drop_na()
