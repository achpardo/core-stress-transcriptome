# Purpose: GO enrichment on core stress genes using the elim algorithm
# Author: Anna Pardo
# Date initiated: Oct. 13, 2023

###### 1. Load modules ########
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)

##### 2. Load data #####
# load gene lists
core_up_setops <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_upgenes_fromsets.txt",header = FALSE)
core_down_setops <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_downgenes_fromsets.txt",header = FALSE)
core_all_setops <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_all_from_setops.txt",header = FALSE)

core_down_total <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_all_downreg.txt",header = FALSE)
core_up_total <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_all_upreg.txt",header = FALSE)
core_all_total <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/union_core_stress_genes.txt",header = FALSE)

core_up_rf <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_fromRF_upreg.txt",header = FALSE)
core_down_rf <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_fromRF_downreg.txt",header = FALSE)
core_all_rf <- read.delim("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/core_genes_from_rf.txt",header = FALSE)

# load GO annotation for maize - B73 V5
zmGO <- fromJSON(file="//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/Zmays_B73V5_GO_annotation.json")

# load go.obo, which was downloaded July 18, 2023 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/go.obo")

##### 4. Gene list data wrangling #####
# convert the two gene vectors to a named list
geneslist <- list()
geneslist[[1]] <- core_up_setops$V1
geneslist[[2]] <- core_down_setops$V1
geneslist[[3]] <- core_all_setops$V1
geneslist[[4]] <- core_down_total$V1
geneslist[[5]] <- core_up_total$V1
geneslist[[6]] <- core_all_total$V1
geneslist[[7]] <- core_up_rf$V1
geneslist[[8]] <- core_down_rf$V1
geneslist[[9]] <- core_all_rf$V1
names(geneslist) <- c("Up_SetOps","Down_SetOps","All_SetOps","Down_Total","Up_Total","All_Total","Up_RF","Down_RF","All_RF")

##### 5. GO enrichment function #####
go_enrichment <- function(injson,spgo){
  ### injson is input json file (named list) with genes of interest ###
  ### spgo is the GO annotation for the species (as named list) ###
  
  allgenes_sp <- names(spgo)
  intgenes <- lapply(injson,function(x){factor(as.integer(allgenes_sp %in% x))})
  for (i in 1:length(intgenes)){
    names(intgenes[[i]]) <- allgenes_sp
  }
  
  all.GOdata.bp <- lapply(intgenes,function(x){
    new("topGOdata",ontology="BP",allGenes=x,annot=annFUN.gene2GO,gene2GO=spgo)})
  
  # initialize test stat - Fisher's exact test
  test.stat <- new("weight01Count",testStatistic=GOFisherTest,name="Fisher test")
  
  all.resultFisher.bp <- lapply(all.GOdata.bp,function(x){getSigGroups(x,test.stat)})
  all.sig.GO.bp <- lapply(all.resultFisher.bp,function(x){data.frame("p-value"=score(x))})
  for (i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]]$GO.term <- row.names(all.sig.GO.bp[[i]])
    all.sig.GO.bp[[i]]$p.adj <- p.adjust(all.sig.GO.bp[[i]]$p.value,method = "fdr")
    all.sig.GO.bp[[i]] <- all.sig.GO.bp[[i]][all.sig.GO.bp[[i]]$p.adj<0.05,]
  }
  all.sig.GO.bp <- all.sig.GO.bp[lapply(all.sig.GO.bp,nrow)>0]
  for (i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]]$regset <- names(all.sig.GO.bp)[i]
    all.sig.GO.bp[[i]]$Ontology <- "BP"
  }
  for(i in 1:length(all.sig.GO.bp)){
    all.sig.GO.bp[[i]] %>% filter(p.adj < 0.05)
  }
  prelimdf <- do.call(rbind,all.sig.GO.bp)
  
  all.sig.GO.res.bp <- lapply(seq_along(all.resultFisher.bp),function(x,y,z){
    GenTable(z[[x]],y[[x]],topNodes = 250)},y=all.resultFisher.bp,z=all.GOdata.bp)
  #add module classification to each df in this list
  for(i in seq_along(all.sig.GO.res.bp)){
    all.sig.GO.res.bp[[i]]$regset <- names(all.resultFisher.bp[i])
  }
  # convert list to df
  sigann_df <- bind_rows(all.sig.GO.res.bp)
  # join to previous df
  dfname <- left_join(prelimdf,sigann_df,by = c("GO.term" = "GO.ID","regset"))
  
  return(dfname)
}

##### 6. Run GO enrichment #####
corego <- go_enrichment(geneslist,zmGO)
# check to see if there are any non-matching GO terms in each dataframe
unique(corego$GO.term[!(corego$GO.term %in% GO$id)])
# the GO terms are all updated!

# now add the descriptions of the GO terms
corego$Descr <- apply(corego,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})

# save final results as .txt file
write.table(corego, file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/GO_enrichment_coregenes_all_weight01.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
