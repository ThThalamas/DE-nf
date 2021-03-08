#if (TPM==TRUE){
#TPM
#  longueur <- read.table("/home/boris/Bureau/projet/projetS2/longueur.txt", row.names = 1, header = F, sep = "\t")
#  for (i in 1:length(rownames(countData))){
#    rate <- log(1+countData[i,])-log(longueur[rownames(countData),])
#    countData[i,] <- exp(rate - log(sum(exp(rate))) + log(1E6))
#  }
#  countData <- round(countData)

#intersect avec les gènes dont on connait la longueur
#  countData <- countData[which(rownames(countData) %in% intersect(rownames(longueur),rownames(countData))),]
#}
#prendre que les femelles, topgo goseq 
#cluster pedago
# pipeline --genomessparsd 4 

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