library(org.Hs.eg.db)
library(Matrix)
library(patchwork)
library(plyr)
library(Seurat)
library(stringr)
library(tidyselect)
library(igraph)
library(pROC)
library(dplyr)


setwd("~/BfBio/")


files<-list.files("raw_omics/")
studies<-c()
for (f in files){
  study<-unlist(strsplit(f , "_") )[1]
  studies<-c(studies , study)
}
studies<-unique(studies)
outdir_data<-paste0("datasets/")


if (!file.exists(outdir_data)){
  dir.create(outdir_data)
  print(paste0(outdir_data , " is created"))
}
outdir_metadata<-paste0("datasets/")

if (!file.exists(outdir_metadata)){
  dir.create(outdir_metadata)
  print(paste0(outdir_metadata , " is created"))
}

for (study in studies){
  idxs<-which(grepl(study , files))
  counts<-data.frame()
  metadata<-data.frame()
  samples<-c()
  for (idx in idxs){
    file_name<-files[idx]
    s<-unlist(strsplit(file_name , "_"))
    if(length(s) == 3 ){
      metadata_file_name<-paste0(s[1] , "_" , s[2], " metadata.csv")  
    }else{
      metadata_file_name<-paste0(s[1] ,  " metadata.csv")  
    }
    
    new_counts<-read.csv(paste0("raw_omics/" , file_name )) 
    new_meta<-read.csv(paste0("raw_omics/" , metadata_file_name ))
    

    
    samples<-c(samples, colnames(new_counts)[colnames(new_counts) != "Feature"])
    
    if (nrow(counts) == 0 ){
      counts<-new_counts
      metadata<-new_meta
    }else{
      counts<-merge(counts , new_counts , by.x = 'Feature' , by.y = 'Feature')
      metadata<-rbind(metadata, new_meta)
    }
  }
  # Map ENTREZIDs to SYMBOLS
  counts$SYMBOL <- mapIds(
    org.Hs.eg.db,
    keys = as.character(counts$Feature),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"  # Take the first match if there are multiple SYMBOLs
  )
  ### check how many symbols are not mapped?
  # Handle duplicates and NAs
  print(paste(sum(is.na(counts$SYMBOL))/nrow(counts) , " ENTREZ are not mapped :-(  for study " , study ) )
  counts <- counts[!duplicated(counts$SYMBOL) & !is.na(counts$SYMBOL), ]
  print(paste("study has the following dimensions " , dim(counts )))
  counts$Feature<-counts$SYMBOL
  if(sum(duplicated(samples)) > 0 ){
    cols = c( "Feature" , samples)
    idx<-which(duplicated(cols))
    counts<-counts[ , -idx]
    print(paste("For study " , study , " removed " , length(idx) , " duplicated samples " ) )
  }
 
  counts<-counts[, colnames(counts) != "SYMBOL"]
  if ("X" %in% colnames(counts)){
    counts<-counts[ , -which(colnames(counts) == "X")]
  }
  write.csv(file = paste0(outdir_data , study, "_data.csv") , counts , row.names = FALSE)
  write.csv(file = paste0(outdir_metadata , study, "_metadata.csv") , metadata , row.names = FALSE)
  
}
