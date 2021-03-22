library(gridExtra)
library(grid)

DE <- DEA_analysis(matrix,metadata)

DE_male_Conventional_vs_Axenic <- DEA_analysis(matrix,metadata)
DE_male_Gnotobiotic_vs_Axenic <- DEA_analysis(matrix,metadata)
DE_male_Gnotobiotic_vs_Conventional <- DEA_analysis(matrix,metadata)
DE_female <- DEA_analysis(matrix,metadata)
DE_female_Conventional_vs_Axenic <- DEA_analysis(matrix,metadata)
DE_female_Gnotobiotic_vs_Axenic <- DEA_analysis(matrix,metadata)
DE_female_Gnotobiotic_vs_Conventional <- DEA_analysis(matrix,metadata)


result <- DEA_result(dataset)
result_male_Conventional_vs_Axenic <- DEA_result(dataset)
result_male_Gnotobiotic_vs_Axenic <- DEA_result(dataset)
result_male_Gnotobiotic_vs_Conventional <- DEA_result(dataset)
result_female <- DEA_result(dataset)
result_female_Conventional_vs_Axenic <- DEA_result(dataset)
result_female_Gnotobiotic_vs_Axenic <- DEA_result(dataset)
result_female_Gnotobiotic_vs_Conventional <- DEA_result(dataset)


DE_male <- DEA_analysis(matrix,metadata)
result_male <- DEA_result(DE_male)

save()

# Rapport
pdf("Cross_Result.pdf", height = 12, width = 12)

#grid.text("Intersect : result male and female", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
#grid.table(intersect(rownames(result_male),rownames(result_female)))
#grid.newpage()

grid.text("All", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result)
rld <- rlog(DE, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Female", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_female)
rld <- rlog(DE_female, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Female : Conventional vs Axenic", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_female_Conventional_vs_Axenic)
rld <- rlog(DE_female_Conventional_vs_Axenic, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Female : Gnotobiotic vs Axenic", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_female_Gnotobiotic_vs_Axenic)
rld <- rlog(DE_female_Gnotobiotic_vs_Axenic, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Female : Gnotobiotic vs Conventional", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_female_Gnotobiotic_vs_Conventional)
rld <- rlog(DE_female_Gnotobiotic_vs_Conventional, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()


grid.text("Male", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_male)
rld <- rlog(DE_male, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Male : Conventional vs Axenic", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_male_Conventional_vs_Axenic)
rld <- rlog(DE_male_Conventional_vs_Axenic, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Male : Gnotobiotic vs Axenic", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_male_Gnotobiotic_vs_Axenic)
rld <- rlog(DE_male_Gnotobiotic_vs_Axenic, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

grid.text("Male : Gnotobiotic vs Conventional", x = 0.05, hjust = 0.05, vjust = -3, gp = gpar(fontsize = 50))
grid.table(result_male_Gnotobiotic_vs_Conventional)
rld <- rlog(DE_male_Gnotobiotic_vs_Conventional, blind = T)
plotPCA(rld, intgroup = "Modality") + theme_bw() + ggtitle("Rlog transformed counts")
grid.newpage()

dev.off()
