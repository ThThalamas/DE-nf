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

if (TPM==TRUE){
  #TPM
  longueur <- read.table("/home/boris/Bureau/projet/projetS2/longueur.txt", row.names = 1, header = F, sep = "\t")
  for (i in 1:length(rownames(countData))){
    rate <- log(1+countData[i,])-log(longueur[rownames(countData),])
    countData[i,] <- exp(rate - log(sum(exp(rate))) + log(1E6))
  }
  countData <- round(countData)
  
  #intersect avec les gènes dont on connait la longueur
  countData <- countData[which(rownames(countData) %in% intersect(rownames(longueur),rownames(countData))),]
}

#prendre que les femelles, topgo goseq 
#cluster pedago
# pipeline --genomessparsd 4 
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
  colData <- colData[which(colData$`Sex`=='Male'),]
  
  countData <- countData[,which(colnames(countData) %in% intersect(colData$`Name of sample`, colnames(countData)))]
  
  #DESeq2 analysis
  dataset <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~Modality)# + Sex
  dataset <- estimateSizeFactors(dataset) #+ dataset@colData$sizeFactor
  
  dataset <- DESeq(dataset)

  return(dataset)
}
matrix <- "/home/boris/Bureau/projet/projetS2/script/finale.txt"
metadata <- "/home/boris/Bureau/projet/projetS2/Metadata.xls"
dataset <- DEA_analysis(matrix, metadata)

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


# PCA
rld <- rlog(dataset, blind = T)
P <- plotPCA(rld, intgroup = "Modality")
P + theme_bw() + ggtitle("Rlog transformed counts")

P <- plotPCA(rld, intgroup = "Sex")
P + theme_bw() + ggtitle("Rlog transformed counts")






############################################################################################################################################### ##
## -- Analyses & Visualisation -- ## 
####################################
# Main result :
View(topResults)
dim(topResults)

# Figure time :
#aheatmap
top40Genes <- head(assay(rld)[rownames(topResults),], 40)
aheatmap(top40Genes, Colv = TRUE, hclustfun = "average", scale = "row")
aheatmap(topGenes, Colv = TRUE, hclustfun = "average", scale = "row") 

# Histograme of p_values
hist(result$pvalue, col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")
hist(result$pvalue[result$baseMean > 1], col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")

# PlotMA
plotMA(topGenes, main = "MA plot: group1 vs group2 (padj < 0.05)", ylim = c(-4,4))
abline(h = c(-1,1), col = "blue", lty = 2)
mtext(c(paste("-2 fold"), paste("+ 2 fold")), side = 4, at = c(-1, 1), cex = 0.8, line = 0.5, col = "blue")

# Dendrogramme
rlog.norm.counts <- assay(rld)
distance.m_rlog  <- as.dist(1 - cor(rlog.norm.counts , method = "pearson" ))
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),main = "rlog  transformed  read  counts\ndistance: Pearson  correlation")


# Plot count
topGene <- rownames(result)[which.min(result$padj)]
els <- plotCounts(dataset, gene = topGene, intgroup="group", returnData=TRUE)
els$sex = c("F","F","F","F","F","F","F","M","M","M","M","M","M","M",
        "F","F","F","F","F","F","F","M","M","M","M","M","M","M",
        "F","F","F","F","M","M","M","M",
        "F","F","F","F","F","F","F","M","M","M","M","M","M","M")

ggplot(els,aes(x = group, y = count)) + ggtitle(topGene) + 
  geom_point(aes(color = sex))

# Other
plotDispEsts(dataset)
