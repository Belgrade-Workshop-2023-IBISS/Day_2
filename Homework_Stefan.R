#model2 analysed in Seurat4
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(ggplot2))

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("dtm2451/dittoSeq")

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
seurat_ex <- readRDS("E:/Bioinformatics/model2.rds")
seurat_ex <- UpdateSeuratObject(seurat_ex)
DimPlot(seurat_ex, label = T) + NoLegend()
#DimPlot(seurat_ex) + scale_color_manual(values = colorblind_vector(2))
gene.sets <- getGeneSets(library = "C5", gene.sets= c("GOBP_MONOCYTE_CHEMOTAXIS", "GOBP_MONOCYTE_DIFFERENTIATION", "GOBP_MONOCYTE_EXTRAVASATION"))
gene.sets                         
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
ridgeEnrichment(ES2, gene.set = "GOBP_MONOCYTE_CHEMOTAXIS", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "GOBP_MONOCYTE_DIFFERENTIATION", add.rug = TRUE)
ridgeEnrichment(ES2, gene.set = "GOBP_MONOCYTE_EXTRAVASATION", add.rug = TRUE)
statistic <- getSignificance(ES2, group = "cluster", fit = "ANOVA") 
statistic
