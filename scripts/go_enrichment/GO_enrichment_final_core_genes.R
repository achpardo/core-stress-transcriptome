# Purpose: GO enrichment on final core stress genes using the elim and weight01 algorithms
# Author: Anna Pardo
# Date initiated: Mar. 25, 2024

###### 1. Load modules ########
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)

##### 2. Load data #####
# load go.obo, which was downloaded Mar. 25, 2024 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/go.obo")

# load GO annotation for maize - B73 V5
zmGO <- fromJSON(file="//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/Zmays_B73V5_GO_annotation.json")

# load lists of core genes
cgall <- fromJSON(file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/allsamp_core_genes_25-Mar-2024.json")
cgpsyn <- fromJSON(file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/psyn_core_genes_25-Mar-2024.json")

##### 4. GO enrichment function #####
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

##### 5. Run GO enrichment for cgall #####
allcorego <- go_enrichment(cgall,zmGO)
# check to see if there are any non-matching GO terms in each dataframe
unique(allcorego$GO.term[!(allcorego$GO.term %in% GO$id)])
# the GO terms are all updated!

# now add the descriptions of the GO terms
allcorego$Descr <- apply(allcorego,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})

# save final results as .txt file
write.table(allcorego, file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/GO_enrichment_finalcoregenes_all_weight01.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)

###### 6. Run GO enrichment for cgpsyn step-by-step ####
allgenes_sp <- names(zmGO)
intgenes <- lapply(cgpsyn,function(x){factor(as.integer(allgenes_sp %in% x))})
for (i in 1:length(intgenes)){
  names(intgenes[[i]]) <- allgenes_sp
}

all.GOdata.bp <- lapply(intgenes,function(x){
  new("topGOdata",ontology="BP",allGenes=x,annot=annFUN.gene2GO,gene2GO=zmGO)})

# initialize test stat - Fisher's exact test
test.stat <- new("weight01Count",testStatistic=GOFisherTest,name="Fisher test")

all.resultFisher.bp <- lapply(all.GOdata.bp,function(x){getSigGroups(x,test.stat)})
all.sig.GO.bp <- lapply(all.resultFisher.bp,function(x){as.data.frame(score(x))})
for (i in 1:length(all.sig.GO.bp)){
  all.sig.GO.bp[[i]]$GO.term <- row.names(all.sig.GO.bp[[i]])
  all.sig.GO.bp[[i]]$p.adj <- p.adjust(all.sig.GO.bp[[i]]$`score(x)`,method = "fdr")
  all.sig.GO.bp[[i]] <- all.sig.GO.bp[[i]][all.sig.GO.bp[[i]]$p.adj<0.05,]
}
all.sig.GO.bp <- all.sig.GO.bp[lapply(all.sig.GO.bp,nrow)>0]
for (i in 1:length(all.sig.GO.bp)){
  all.sig.GO.bp[[i]]$regset <- names(all.sig.GO.bp)[i]
  #all.sig.GO.bp[[i]]$Ontology <- "BP"
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
