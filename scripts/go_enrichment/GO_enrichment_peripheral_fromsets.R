# Purpose: GO enrichment on putative peripheral stress responsive genes (unique to each stressor) from set operations.
# Author: Anna Pardo
# Date initiated: July 30, 2023

###### 1. Load modules ########
library(topGO)
library(rjson)
library(dplyr)
library(ontologyIndex)

##### 2. Load data #####
# load upregulated and downregulated gene sets
upgenes <- fromJSON(file = "//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/peripheral_response_upgenes_fromsets.json")
downgenes <- fromJSON(file = "//wsl.localhost/Ubuntu/home/leviathan22/core-stress-transcriptome/data/peripheral_response_downgenes_fromsets.json")

# load GO annotation
zmgo <- fromJSON(file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/Zmays_B73V5_GO_annotation.json")

# load go.obo, which was downloaded July 18, 2023 from http://geneontology.org/docs/download-ontology/#go_obo_and_owl
GO <- get_ontology("//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/go.obo")

##### 3. GO enrichment function #####
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
  test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test")
  
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

##### 4. Run GO enrichment #####
upgo <- go_enrichment(upgenes,zmgo)
downgo <- go_enrichment(downgenes,zmgo)

# check to see if there are any non-matching GO terms in each dataframe
unique(upgo$GO.term[!(upgo$GO.term %in% GO$id)])
unique(downgo$GO.term[!(downgo$GO.term %in% GO$id)])
downgo$GO.term <- gsub("GO:0006464","GO:0036211",downgo$GO.term)

# now add the descriptions of the GO terms
upgo$Descr <- apply(upgo,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})
downgo$Descr <- apply(downgo,MARGIN=1,function(x){
  get_term_property(GO,"name",x[2])})

# save final results as .txt file
write.table(upgo, file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/GO_enrichment_peripheral_upgenes_fromsets.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)
write.table(downgo, file = "//wsl$/Ubuntu/home/leviathan22/core-stress-transcriptome/data/GO_enrichment_peripheral_downgenes_fromsets.txt",
            sep = "\t",row.names = FALSE,col.names = TRUE)