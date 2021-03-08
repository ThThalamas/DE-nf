# TODOUX
# prendre que les femelles
# topgo goseq 
# pipeline --genomessparsd 4 

############################################################################################################################################### ##
## -- Analyses & Visualisation -- ## 
####################################
# Main result :
View(topResults)
dim(topResults)

# Figure time :
# Histograme of p_values
hist(result$pvalue, col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")
hist(result$pvalue[result$baseMean > 1], col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies  of p-values")

# Other
plotDispEsts(dataset)


#TPM
if (TPM==TRUE){
  longueur <- read.table("/home/boris/Bureau/projet/projetS2/longueur.txt", row.names = 1, header = F, sep = "\t")
  for (i in 1:length(rownames(countData))){
    rate <- log(1+countData[i,])-log(longueur[rownames(countData),])
    countData[i,] <- exp(rate - log(sum(exp(rate))) + log(1E6))
  }
  countData <- round(countData)
  
  #intersect avec les gènes dont on connait la longueur
  countData <- countData[which(rownames(countData) %in% intersect(rownames(longueur),rownames(countData))),]
}