# Homework week 2

library(escape)
library(Seurat)
library(ggplot2)
library(dittoSeq)

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

seurat_ex <- readRDS("D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week2/homework2/pbmc3k_final.rds")
seurat_ex <- UpdateSeuratObject(seurat_ex)
DimPlot(seurat_ex, label = T) + NoLegend()

library(escape)

gene.sets <- getGeneSets(library = "C5", gene.sets =  c("GOBP_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION", "GOBP_MONONUCLEAR_CELL_MIGRATION", "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION"))

ES <- enrichIt(obj = seurat_ex, gene.sets = gene.sets, groups = 1000)
seurat_ex <- AddMetaData(seurat_ex, ES)

seurat_ex@meta.data$active.idents <- seurat_ex@active.ident

dittoHeatmap(seurat_ex, genes = NULL, metas = names(ES),
             heatmap.colors = rev(colorblind_vector(50)),
             annot.by = c("active.idents"),
             cluster_cols = TRUE,
             fontsize = 7)


ES2 <- data.frame(seurat_ex[[]], Idents(seurat_ex))
colnames(ES2)[ncol(ES2)] <- "cluster"

ridgeEnrichment(ES2, gene.set = "GOBP_REGULATION_OF_MONONUCLEAR_CELL_MIGRATION", group = "cluster", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "GOBP_MONONUCLEAR_CELL_MIGRATION", group = "cluster", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "GOBP_MONONUCLEAR_CELL_DIFFERENTIATION", group = "cluster", add.rug = TRUE)

statistics <- getSignificance(ES2, group = "cluster", fit = "ANOVA") 




