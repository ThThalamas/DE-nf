#!/usr/bin/Rscript
###############################################################################################################################################
library("DESeq2")
library ("edgeR")
library("pheatmap")
library("RColorBrewer")
library("ggbeeswarm")
library("genefilter")
library(biomaRt)
library(stringr)
library(ggplot2)
library(NMF)
library(tidyverse)

remplace <- function(name){name <- gsub("-", ".", name)} 
split <- function(name){name <- str_split(name,"_")[[1]][1]}


###############################################################################################################################################
## -- DESeq2 annalyse -- ##
###########################
# DEA
DEA_analysis <- function(matrix, metadata){
  #load
  countData <- read.table(matrix, row.names = 1, header = T, sep = "\t")
  countData <- head(countData,-5)
  countData <- countData[1:61]
  
  #cutoff
  cutoff <- as.integer(cpm(10, mean(colSums(countData))))
  keep <- rowSums(cpm(countData)>cutoff) >= 2
  countData <- countData[keep,]
  
  #rename
  colnames(countData) <- apply(as.matrix(colnames(countData)), 1, split)
  countData <- countData[,order(colnames(countData))]
  
  #process metadata
  colData <- readxl::read_xls(path = metadata)
  colData[1] <- apply(colData[1], 1, remplace)
  colData <- colData[order(colData[,1]),]
  
  #intersect with data laready processed
  colData <- colData[which(colData$`Name of sample` %in% intersect(colData$`Name of sample`, colnames(countData))),]
  #colData <- colData[which(colData$`Sex`=='Female'),]
  
  countData <- countData[,which(colnames(countData) %in% intersect(colData$`Name of sample`, colnames(countData)))]
  
  #DESeq2 analysis
  dataset <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~Modality)# + Sex
  dataset <- estimateSizeFactors(dataset) #+ dataset@colData$sizeFactor
  
  dataset <- DESeq(dataset)

  return(dataset)
}
#matrix <- "/home/boris/Bureau/projet/DE-nf/result/finale.txt"
#metadata <- "/home/boris/Bureau/projet/DE-nf/data/Metadata.xls"

matrix = as.character(commandArgs(TRUE)[1])
metadata = as.character(commandArgs(TRUE)[2])
dataset <- DEA_analysis(matrix, metadata)
rld <- rlog(dataset, blind = T)



# PCA
pdf("/home/boris/Bureau/projet/DE-nf/result/Result_PCA.pdf")
P <- plotPCA(rld, intgroup = "Modality")
P + theme_bw() + ggtitle("Rlog transformed counts")

P <- plotPCA(rld, intgroup = "Sex")
P + theme_bw() + ggtitle("Rlog transformed counts")
dev.off()



# Top Genes Results
DEA_result <- function(dataset, group, reference, group2){
  if (group != "Sex"){result <- results(dataset, c(group, reference, group2), pAdjustMethod = "BH")}
  else{result <- results(dataset, pAdjustMethod = "BH")}
  result <- result[complete.cases(result),]
  result <- result[order(result$padj),]
  
  #Tri des gènes les plus différentiellement exprimé
  topResults <- as.matrix(subset(result , padj < 0.05))
  
  return(topResults)
}
#By Modality
modality_result <- DEA_result(dataset, "Modality" , "Axenic" , "Conventional") #"Gnotobiotic Microbacterium oxydans", "Gnotobiotic Carnobacterium maltaromatycum", "Gnotobiotic Oerskovia"
topGenes_modality <- assay(rld)[rownames(modality_result),]

#By Sex
sex_result <- DEA_result(dataset, "Sex", "Female", "Male")
topGenes_sex <- assay(rld)[rownames(sex_result),]


