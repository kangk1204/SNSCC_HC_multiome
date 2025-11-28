library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# load RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")

snscc_filtered$CellType_sample <- paste0(snscc_filtered$CellType,"(",snscc_filtered$sample,")")
table(snscc_filtered$CellType_sample)

snscc_filtered$CellType_dataset <- paste0(snscc_filtered$CellType,"(",snscc_filtered$dataset,")")
table(snscc_filtered$CellType_dataset)

color.snscc.filtered.CellType <- c("#d42a34", "#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")
color.sample <- c("#F7C475","#7D5008")

Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")

DimPlot(snscc_filtered, label = F, reduction = "umap.wnn", cols = color.sample, pt.size = 0.01, group.by = "sample") -> Figure1D
Figure1D
ggsave(filename = "./Figures/snscc_filtered_res_08_sample_wnnUMAP_5*4.pdf", plot = Figure1D, width = 5, height = 4, units = "in")

DimPlot(snscc_filtered, label = F, reduction = "umap.wnn", cols = color.snscc.filtered.CellType, pt.size = 0.01) -> Figure1E-a
Figure1E-a
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_wnnUMAP_5.7*4.pdf", plot = Figure1E-a, width = 5.7, height = 4, units = "in")

DimPlot(snscc_filtered, label = F, reduction = "umap.rna", cols = color.snscc.filtered.CellType, pt.size = 0.01) -> Figure1E-b
Figure1E-b
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_rnaUMAP_5.7*4.pdf", plot = Figure1E-b, width = 5.7, height = 4, units = "in")

DimPlot(snscc_filtered, label = F, reduction = "umap.atac", cols = color.snscc.filtered.CellType, pt.size = 0.01) -> Figure1E-c
Figure1E-c
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_atacUMAP_5.7*4.pdf", plot = Figure1E-c, width = 5.7, height = 4, units = "in")



# dittoSeq
library(dittoSeq)

dittoBarPlot(snscc_filtered, var = "CellType", group.by = "sample", color.panel = color.snscc.filtered.CellType, var.labels.reorder = c(2,3,1,4,6,5), x.labels.rotate = F) -> Figure1G
Figure1G
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_proportion_4*3.pdf", plot = Figure1G, width = 3, height = 4, units = "in")



# scProportionTest
library("scProportionTest")

prop_test <- sc_utils(snscc_filtered)

prop_test <- permutation_test(
  prop_test, cluster_identity = "CellType",
  sample_1 = "HC", sample_2 = "Tumor",
  sample_identity = "sample")

permutation_plot(prop_test, FDR_threshold = 0.05, log2FD_threshold = log2(2), order_clusters = TRUE) -> FigureS1G
FigureS1G
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_scProportionTest_6*3.pdf", plot = FigureS1G, width = 6, height = 3, units = "in")



# Marker gene expression & track
# FeaturePlot
library(RColorBrewer)
library(scCustomize)

FeaturePlot_scCustom(snscc_filtered, features = c("CDH1"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_CDH1_gene
Figure1F_CDH1_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_CDH1_4.5*4.pdf", plot = Figure1F_CDH1_gene, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("PDGFRB"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_PDGFRB_gene
Figure1F_PDGFRB_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_PDGFRB_4.5*4.pdf", plot = Figure1F_PDGFRB_gene, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("VWF"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_VWF_gene
Figure1F_VWF_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_VWF_4.5*4.pdf", plot = Figure1F_VWF_gene, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("ITGAX"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_ITGAX_gene
Figure1F_ITGAX_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_ITGAX_4.5*4.pdf", plot = Figure1F_ITGAX_gene, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("CD2"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_CD2_gene
Figure1F_CD2_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_CD2_4.5*4.pdf", plot = Figure1F_CD2_gene, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("FCRL5"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> Figure1F_FCRL5_gene
Figure1F_FCRL5_gene
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_FCRL5_4.5*4.pdf", plot = Figure1F_FCRL5_gene, height = 4, width = 4.5, units = "in")

# Track
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

snscc_filtered_regionstats <- RegionStats(snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

Idents(snscc_filtered_regionstats) <- "CellType"

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CDH1")
CoveragePlot(object = snscc_filtered_regionstats, region = "CDH1", extend.upstream = 6000, extend.downstream = -90000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_CDH1
Figure1F_CDH1
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_CDH1_4*6.pdf", plot = Figure1F_CDH1, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "PDGFRB")
CoveragePlot(object = snscc_filtered_regionstats, region = "PDGFRB", extend.upstream = -36000, extend.downstream = 8000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_PDGFRB
Figure1F_PDGFRB
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_PDGFRB_4*6.pdf", plot = Figure1F_PDGFRB, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "VWF")
CoveragePlot(object = snscc_filtered_regionstats, region = "VWF", extend.upstream = -150000, extend.downstream = 20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_VWF
Figure1F_VWF
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_VWF_4*6.pdf", plot = Figure1F_VWF, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "ITGAX")
CoveragePlot(object = snscc_filtered_regionstats, region = "ITGAX", extend.upstream = 6000, extend.downstream = -20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_ITGAX
Figure1F_ITGAX
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_ITGAX_4*6.pdf", plot = Figure1F_ITGAX, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CD2")
CoveragePlot(object = snscc_filtered_regionstats, region = "CD2", extend.upstream = 18000, extend.downstream = -10000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_CD2
Figure1F_CD2
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_CD2_4*6.pdf", plot = Figure1F_CD2, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "FCRL5")
CoveragePlot(object = snscc_filtered_regionstats, region = "FCRL5", extend.upstream = 10000, extend.downstream = -20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> Figure1F_FCRL5
Figure1F_FCRL5
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_Track_FCRL5_4*6.pdf", plot = Figure1F_FCRL5, height = 6, width = 4, units = "in")



# DEGs (CellType_sample; FindAllMarkers)
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

DefaultAssay(snscc_filtered) <- "RNA"

snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05 & snscc_filtered.all.genes$avg_log2FC > 1, ]
write.table(snscc_filtered.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_DEGs_log2FC_1_HC_Tumor.txt", sep = '\t')

snscc_filtered.da.genes.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.genes$gene, assays = "RNA")
pheatmap::pheatmap(snscc_filtered.da.genes.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> Figure2A
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_sample_DEGs_Log2FC_1_heatmap_2*3.5.pdf", plot = Figure2A, width = 2, height = 3.5, units = "in")

gene_counts <- table(snscc_filtered.all.da.genes$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType_sample", ylab = "Gene number")



# Gene expression (CellType_sample)
# DotPlot
DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
DotPlot(snscc_filtered, features = c("ADM","ENO1","LDHA","NDRG1","PKM","SLC2A1","CDH3","DSC3","DSG3","KRT6A","KRT17","SFN","FAP","POSTN","COL1A1","COL1A2","COL6A1","TGFB2","TGFBR1","ITGAV","ITGB1","ACVRL1","ANGPT2","FLT4","KDR","NOTCH1","NOTCH4","ACP5","APOE","APOC1","ADAM8","AIF1","NFKBID","ABCA1","ICAM1","NAMPT","TNFRSF14","TNFRSF1B","TNFRSF18","ICOS","IL2RB","JAK3"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip() -> Figure2B
Figure2B
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_5.5*10.5.pdf", plot = Figure2B, width = 5.5, height = 10.5, units = "in")

DotPlot(snscc_filtered, features = c("CCL28","DMBT1","EGF","STATH"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + theme(axis.title.x = element_blank()) -> FigureS2C_HC_epi
FigureS2C_HC_epi
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_HC_Epi_4.5*4.5.pdf", plot = FigureS2C_HC_epi, width = 4.5, height = 4.5, units = "in")

DotPlot(snscc_filtered, features = c("ACTA2","FHL1","MYH11","MYOCD"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + theme(axis.title.x = element_blank()) -> FigureS2C_HC_Fibro-b
FigureS2C_HC_Fibro
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_HC_Fibro_4.5*4.5.pdf", plot = FigureS2C_HC_Fibro, width = 4.5, height = 4.5, units = "in")

# FeaturePlot
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("KRT6A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure2D_KRT6A
Figure2D_KRT6A
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_KRT6A_4*9.pdf", plot = Figure2D_KRT6A, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NDRG1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure2D_NDRG1
Figure2D_NDRG1
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_NDRG1_4*9.pdf", plot = Figure2D_NDRG1, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FAP"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure2D_FAP
Figure2D_FAP
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_FAP_4*9.pdf", plot = Figure2D_FAP, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ANGPT2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure2D_ANGPT2
Figure2D_ANGPT2
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_ANGPT2_4*9.pdf", plot = Figure2D_ANGPT2, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("IL2RB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure2D_IL2RB
Figure2D_IL2RB
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_IL2RB_4*9.pdf", plot = Figure2D_IL2RB, height = 4, width = 9, units = "in")



# Progeny
# We create a data frame with the specification of the cells that belong to 
library(progeny)
library(magrittr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

## each cluster to match with the Progeny scores.
DefaultAssay(snscc_filtered) <- "RNA"

Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

snscc_filtered$sample_CellType <- paste0(snscc_filtered$sample,"_",snscc_filtered$CellType)

Idents(snscc_filtered) <- "sample_CellType"

snscc_filtered <- RenameIdents(snscc_filtered, "HC_Epithelial cells"="A","HC_Fibroblasts"="B","HC_Endothelial cells"="C","HC_Myeloid cells"="D","HC_T cells"="E","HC_Plasma cells"="F","Tumor_Epithelial cells"="G","Tumor_Fibroblasts"="H","Tumor_Endothelial cells"="I","Tumor_Myeloid cells"="J","Tumor_T cells"="K","Tumor_Plasma cells"="L")
snscc_filtered$progeny_cell_order <- Idents(snscc_filtered)
Idents(snscc_filtered) <- "progeny_cell_order"

CellsClusters <- data.frame(Cell = names(Idents(snscc_filtered)), CellType = as.character(Idents(snscc_filtered)), stringsAsFactors = FALSE)
DimPlot(snscc_filtered, label = TRUE, pt.size = 0.01, reduction = "umap.wnn") + NoLegend()

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 
snscc_filtered_progeny <- progeny(snscc_filtered, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
snscc_filtered_progeny <- Seurat::ScaleData(snscc_filtered_progeny, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(snscc_filtered_progeny, slot = "scale.data", assay = "progeny"))) %>% rownames_to_column("Cell") %>% gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(11,"RdBu")))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 90,
                        treeheight_col = 20,  border_color = "black", treeheight_row = 20, cutree_rows = 4, cluster_cols = F) -> Figure2C
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_sample_PROGENy_1.5*5.pdf", plot = Figure2C, width = 1.5, height = 5, units = "in")


##########################################################





























































## AddModuleScore
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scCustomize)

DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
AddModuleScore_gene_list <- read.csv("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Data/AddModuleScore_2024.csv", header = T, sep = ",")

Module <- list(AddModuleScore_gene_list$MsigDB_hypoxia)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "A")
FeaturePlot_scCustom(snscc_filtered, features = "A1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "OrRd")))

Module <- list(AddModuleScore_gene_list$bulk_keratinocyte_diff)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "bulk_keratinocyte_diff")
VlnPlot(snscc_filtered, features = "bulk_keratinocyte_diff1", pt.size = 0, cols = color.snscc.filtered.CellType, idents = c("Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/Modulescore_bulk_keratinocyte_diff_vlnplot_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")

#gene_list <- unique(AddModuleScore_gene_list$bulk_chromatin_remodeling)
#DotPlot(snscc_filtered, features = gene_list, cols = c("#CCCCCC","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev)

Module <- list(AddModuleScore_gene_list$bulk_PID_HIF1)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "bulk_PID_HIF1")
VlnPlot(snscc_filtered, features = "bulk_PID_HIF11", pt.size = 0, cols = color.snscc.filtered.CellType, idents = c("Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/Modulescore_bulk_PID_HIF_vlnplot_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")

Module <- list(AddModuleScore_gene_list$bulk_Proteoglycans_in_cancer)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "bulk_proteo")
VlnPlot(snscc_filtered, features = "bulk_proteo1", pt.size = 0, cols = color.snscc.filtered.CellType, idents = c("Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/Modulescore_bulk_proteo_vlnplot_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")

Module <- list(AddModuleScore_gene_list$bulk_VEGFA)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "VEGF")
VlnPlot(snscc_filtered, features = "VEGF1", pt.size = 0, cols = color.snscc.filtered.CellType, idents = c("Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/Modulescore_bulk_VEGF_vlnplot_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")


AddModuleScore_gene_list <- read.csv("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Data/AddModuleScore2.csv", header = T, sep = ",")
Module <- list(AddModuleScore_gene_list$Basal)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "Basal")
FeaturePlot_scCustom(snscc_filtered, features = "Basal1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "Reds")))

Module <- list(AddModuleScore_gene_list$Classical)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "Classical")
FeaturePlot_scCustom(snscc_filtered, features = "Classical1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "Reds")))

AddModuleScore_gene_list <- read.csv("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Data/AddModuleScore_2024.csv", header = T, sep = ",")
Module <- list(AddModuleScore_gene_list$bulk_down_DEGs_top100)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "Down")
FeaturePlot_scCustom(snscc_filtered, features = "Down1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "PuBu"))) -> p1
Module <- list(AddModuleScore_gene_list$bulk_up_DEGs_top100)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "Up")
FeaturePlot_scCustom(snscc_filtered, features = "Up1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "OrRd"))) -> p2
library(ggpubr)
ggarrange(p1,p2,nrow=2) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_bulk_DEGs_Top_100_module_score_9.5*8.pdf", plot = figure, height = 8, width = 9.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = "Down1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "PuBu")))
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_bulk_down_DEGs_Top_100_module_score_9.5*4.pdf", plot = figure, height = 4, width = 9.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = "Up1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample", na_cutoff = NULL, colors_use = (brewer.pal(n=9, name = "OrRd")))
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_bulk_up_DEGs_Top_100_module_score_9.5*4.pdf", plot = figure, height = 4, width = 9.5, units = "in")

AddModuleScore_gene_list <- c("ABCB1","ABCC1","ABCC2","ABCC3","ABCC5","ABCG2","AHR","AP1S1","APC","AR","ARNT","ATM","BAX","BCL2","BCL2L1","BLMH","BRCA1","BRCA2","CCND1","CCNE1","CDK2","CDK4","CDKN1A","CDKN1B","CDKN2A","CDKN2D","CLPTM1L","CYP1A1","CYP1A2","CYP2B6","CYP2C19","CYP2C8","CYP2C9","CYP2D6","CYP2E1","CYP3A4","CYP3A5","DHFR","EGFR","ELK1","EPHX1","ERBB2","ERBB3","ERBB4","ERCC3","ESR1","ESR2","FGF2","FOS","GSK3A","GSTP1","HIF1A","IGF1R","IGF2R","MET","MSH2","MVP","MYC","NFKB1","NFKB2","NFKBIB","NFKBIE","PPARA","PPARD","PPARG","RARA","RARB","RARG","RB1","RELB","RXRA","RXRB","SOD1","SULT1E1","TNFRSF11A","TOP1","TOP2A","TOP2B","TP53","TPMT","UGCG","XPA","XPC","ACTB","B2M","GAPDH","HPRT1","RPLP0")
Module <- list(AddModuleScore_gene_list)
snscc_filtered <- AddModuleScore(object = snscc_filtered, features = Module, name = "Drug_resist")
FeaturePlot_scCustom(snscc_filtered, features = "Drug_resist1", pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#002153","#00366C","#0072B4","#79ABE2","#FFFFBF","#D07C42","#C7522B","#b21704","#611300"), na_cutoff = NULL, split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_Drug_resistance_module_score_10*4.pdf", plot = figure, width = 10, height = 4, units = "in")
FeaturePlot_scCustom(snscc_filtered, features = "Drug_resist1", pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(brewer.pal(n=11,"RdYlBu")), na_cutoff = NULL, split.by = "sample") -> figure
figure

library(ggsignif)
library(ggplot2)
Idents(snscc_filtered) <- "sample"
my_comparisons <- list(c("Normal","Tumor"))
VlnPlot(snscc_filtered, features = "Drug_resist1", pt.size = 0.01, cols = c("#003366","#CC0000")) # for setting ylim
VlnPlot(snscc_filtered, features = "Drug_resist1", pt.size = 0, cols = c("#003366","#CC0000")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_signif(comparison = my_comparisons, y_position = 0.44, textsize = 4, map_signif_level = T) + ylim(-0.2,0.5) + NoLegend() -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_Drug_resistance_module_score_sample_vlnplot_2*3.pdf", plot = figure, width = 2, height = 3, units = "in")

## Gene expression
library(RColorBrewer)
library(scCustomize)
DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Epithelial cells(Tumor)","Fibroblasts(HC)","Fibroblasts(Tumor)","Endothelial cells(HC)","Endothelial cells(Tumor)","Myeloid cells(HC)","Myeloid cells(Tumor)","T cells(HC)","T cells(Tumor)","Plasma cells(HC)","Plasma cells(Tumor)")

VlnPlot(snscc_filtered, features = "HIF1A", pt.size = 0, cols = c("#fe8181","#cb2424","#ffb38a","#ff6700","#FADA5E","#DAA520"), idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)","Fibroblasts(HC)","Fibroblasts(Tumor)","Endothelial cells(HC)","Endothelial cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_vlnplot_3CellType_HIF1A_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")

VlnPlot(snscc_filtered, features = "JUNB", pt.size = 0, cols = c("#fe8181","#cb2424","#ffb38a","#ff6700","#FADA5E","#DAA520"), idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)","Fibroblasts(HC)","Fibroblasts(Tumor)","Endothelial cells(HC)","Endothelial cells(Tumor)")) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() & theme(text = element_text(size = 20)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_vlnplot_3CellType_JUNB_4*3.pdf", plot = figure, height = 3, width = 4, units = "in")

DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered_motif) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NDRG1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NDRG1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("DSG3"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_DSG3_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("CDH3"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_CDH3_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("DSC3"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_DSC3_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("SLC2A1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_SLC2A1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("LDHA"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_LDHA_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("KRT6A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_KRT6A_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("SFN"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_SFN_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ACVRL1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_FAP_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FAP"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_FAP_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("COL1A1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_COL1A1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NOTCH1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NOTCH1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ANGPT2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_ANGPT2_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ACVRL1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_ACVRL1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("APOE"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_APOE_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("APOC1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_APOC1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ADAM8"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_ADAM8_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFKBID"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFKBID_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ABCA1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_ABCA1_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NAMPT"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NAMPT_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("JAK3"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_JAK3_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("IL2RB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_IL2RB_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot(snscc_filtered, features = c("CDH1"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_CDH1_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("ACTA2"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=11,"Reds")))

FeaturePlot(snscc_filtered, features = c("CDH1"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_CDH1_4.6*4.pdf", plot = figure, height = 4, width = 4.6, units = "in")

FeaturePlot(snscc_filtered, features = c("PDGFRB"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_PDGFRB_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot(snscc_filtered, features = c("VWF"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_VWF_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot(snscc_filtered, features = c("ITGAX"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_ITGAX_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot(snscc_filtered, features = c("CD2"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_CD2_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot(snscc_filtered, features = c("FCRL5"), reduction = "umap.wnn", pt.size = 0.01, order = T, cols = c("#DDDDDD","#CC0000"), max.cutoff = "q95") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_FCRL5_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")


FeaturePlot_scCustom(snscc_filtered, features = c("CDH1"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_CDH1_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("PDGFRB"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_PDGFRB_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("VWF"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_VWF_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("ITGAX"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_ITGAX_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("CD2"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_CD2_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, features = c("FCRL5"), reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"YlOrRd"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_scCustom_FCRL5_4.5*4.pdf", plot = figure, height = 4, width = 4.5, units = "in")


# Gene exp check for motif
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("HIF1A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_HIF1A_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FOSL1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_FOSL1_4*9.2.pdf", plot = figure, height = 4, width = 9.2, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("JUNB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_JUNB_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FOSL2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_FOSL2_4*9.2.pdf", plot = figure, height = 4, width = 9.2, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFATC2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(11,"YlOrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFATC2_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("TP63"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(11,"YlOrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_TP63_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("CEBPD"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(11,"YlOrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_CEBPD_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("KLF6"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(11,"YlOrRd")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_KLF6_4*9.pdf", plot = figure, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(brewer.pal(11,"RdYlBu")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFIB_4*9.2.pdf", plot = figure, height = 4, width = 9.2, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIX"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> figure
figure
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIX"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(brewer.pal(11,"RdYlBu")), split.by = "sample") -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFIX_4*9.2.pdf", plot = figure, height = 4, width = 9.2, units = "in")
#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("HIF1A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = colorRampPalette(ArchRPalettes$horizonExtra)(7), split.by = "sample") -> figure
#figure

#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("HIF1A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(hcl.colors(11,"YlOrRd")), split.by = "sample") -> figure
#figure



Idents(snscc_filtered) <- "sample"
DotPlot(snscc_filtered, features = c("DMBT1","MUC7","MYH11","PIP","PRB1","STATH","ADM","DSC3","KRT6A","MKI67","NDRG1"), cols = c("#DDDDDD","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip()  -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_dotplot_6*3.7.pdf", plot = figure, width = 6, height = 3.7, units = "in")

FeaturePlot(snscc_filtered, features = c("TOP2A"), reduction = "umap.wnn", split.by = "dataset", pt.size = 0.05)
VlnPlot(snscc_filtered, features = c("DMBT1","ENPP3","PROM1","STATH","ADM","DSC3","KRT6A","NDRG1"), pt.size = 0, cols = c("#003399","#003399","#003399","#003399","#CC0000","#CC0000","#CC0000","#CC0000"), stack = T, flip = T) -> figure # + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) -> figure # + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) -> figure
figure
VlnPlot(snscc_filtered, features = c("DMBT1","ENPP3","PROM1","STATH","ADM","DSC3","KRT6A","NDRG1"), pt.size = 0, cols = c("#003399","#003399","#003399","#003399","#CC0000","#CC0000","#CC0000","#CC0000"), stack = T, flip = T) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_sample_vlnplot_1.8*6.pdf", plot = figure, width = 1.8, height = 6, units = "in")

Idents(snscc_filtered) <- "CellType_sample"
DotPlot(snscc_filtered, features = c("CDH1","EPCAM","PDGFRB","THY1","CDH5","VWF","AIF1","ITGAX","CD2","CD3D","FCRL5","JCHAIN"), cols = c("#CCCCCC","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_dotplot_6*3.7.pdf", plot = figure, width = 6, height = 3.7, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_dotplot_6*2.pdf", plot = figure, width = 6, height = 2, units = "in")

Idents(snscc_filtered) <- "CellType"
DotPlot(snscc_filtered, features = c("CDH1","KRT8","DCN","PDGFRB","CD34","VWF","ITGAX","SPI1","CD2","CD3D","FCRL5","JCHAIN"), cols = c("#CCCCCC","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip() -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_marker_gene_dotplot_3.8*4.pdf", plot = figure, width = 3.8, height = 4, units = "in")

DotPlot(snscc_filtered, features = c("CDH1","EPCAM","DCN","THY1","CDH5","VWF","AIF1","ITGAX","CD2","CD3D","FCRL5","JCHAIN"), cols = c("#336699","#CC0000"), split.by = "sample") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip() -> figure
figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_splitby_sample_dotplot_6*3.7.pdf", plot = figure, width = 6, height = 3.7, units = "in")

DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
DotPlot(snscc_filtered, features = c("ADM","ENO1","LDHA","NDRG1","PKM","SLC2A1","CDH3","DSC3","DSG3","KRT6A","KRT17","SFN","FAP","POSTN","COL1A1","COL1A2","COL6A1","TGFB2","TGFBR1","ITGAV","ITGB1","ACVRL1","ANGPT2","FLT4","KDR","NOTCH1","NOTCH4","ACP5","APOE","APOC1","ADAM8","AIF1","NFKBID","ABCA1","ICAM1","NAMPT","TNFRSF14","TNFRSF1B","TNFRSF18","ICOS","IL2RB","JAK3"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip() -> figure
figure
#+ geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_5.5*10.5.pdf", plot = figure, width = 5.5, height = 10.5, units = "in")

DotPlot(snscc_filtered, features = c("CCL28","DMBT1","EGF","STATH"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + theme(axis.title.x = element_blank())  -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_HC_Epi_4.5*4.5.pdf", plot = figure, width = 4.5, height = 4.5, units = "in")

DotPlot(snscc_filtered, features = c("ACTA2","FHL1","MYH11","MYOCD"), cols = "RdBu") + theme_linedraw(base_line_size = 0.8, base_rect_size = 0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + theme(axis.title.x = element_blank()) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_GO_gene_dotplot_HC_Fibro_4.5*4.5.pdf", plot = figure, width = 4.5, height = 4.5, units = "in")

FeaturePlot_scCustom(snscc_filtered, reduction = "umap.wnn", features = c("BMP4"), pt.size = 0.01, colors_use = rev(hcl.colors(11,"YlOrRd")), num_columns = 1, split.by = "sample")

FeaturePlot(snscc_filtered, reduction = "umap.wnn", features = "FHL1", split.by = "dataset", pt.size = 0.01)
FeaturePlot(T_cell, reduction = "umap.wnn", features = "BCL2", split.by = "dataset", pt.size = 0.01)
VlnPlot(T_cell,"STAT5B")
VlnPlot(T_cell,"PSME2")

FeaturePlot(Mye, reduction = "umap.wnn", features = "ID2", split.by = "dataset", pt.size = 0.01)
VlnPlot(Mye,"ICAM1")
VlnPlot(Mye,"PSME2")

#c("DDIT4","KRT17","KRT6A","TP63","VEGFA","ADM","ENO1","LDHA","SLC2A1","COL1A1","COL1A2","FN1","MYO1B","EMP1","GJA1","PNP","APOE","CD68","TYROBP","SPP1","BST2","IL32","SMC4","TAP1","CFL1","PSME2","HAVCR2","TIGIT","ICOS","PFN1","RNF213","PPIA","CTSC","S100A9","S100A8")
#DotPlot(snscc_filtered, features = c("DDIT4","KRT17","KRT6A","TP63","VEGFA","ADM","ENO1","LDHA","SLC2A1","COL1A1","COL1A2","FN1","MYO1B","EMP1","GJA1","PNP","APOE","CD68","TYROBP","SPP1","BST2","IL32","TAP1","PSME2"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> figure
#figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_dotplot_8.5*3.pdf", plot = figure, width = 8.5, height = 3, units = "in")
#DotPlot(snscc_filtered, features = c("DDIT4","KRT17","KRT6A","TP63","VEGFA","ADM","ENO1","LDHA","SLC2A1","COL1A1","COL1A2","FN1","MYO1B","EMP1","GJA1","PNP","APOE","CD68","TYROBP","SPP1","BST2","IL32","TAP1","PSME2"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + scale_x_discrete(limits=rev) -> figure
#figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_dotplot_order_change_8.5*3.pdf", plot = figure, width = 8.5, height = 3, units = "in")


DotPlot(snscc_filtered, features = c("ADM","COL17A1","DSG3","ENO1","FERMT1","GPR87","KRT5","KRT6A","NDRG1","S100A10","S100A2","SUB1","CD68","IL32"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> figure
figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_dotplot_8.5*3.pdf", plot = figure, width = 8.5, height = 3, units = "in")

DotPlot(snscc_filtered, features = c("ADM","VEGFA","CALCRL","KDR"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_dotplot_signaling_4.5*3.5.pdf", plot = figure, width = 4.5, height = 3.5, units = "in")

DotPlot(snscc_filtered, features = c("ADM","CALCRL","RAMP1","RAMP2","RAMP3","VEGFA","FLT1","KDR"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> figure
figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_dotplot_signaling2_4.5*3.5.pdf", plot = figure, width = 4.5, height = 3.5, units = "in")
DotPlot(snscc_filtered, features = c("ADM","CALCRL","RAMP1","RAMP2","RAMP3","BCL2","MAPK1","RAF1","PIK3CA","AKT1"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev)

DotPlot(snscc_filtered, features = c("DNMT1","DNMT3A","DNMT3B","UHRF1","PCNA","TET1","TET2","TET3","TDG","SMUG1","KDM5B"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev)
VlnPlot(snscc_filtered, features = c("DNMT1","DNMT3A","DNMT3B","UHRF1","PCNA","TET1","TET2","TET3","TDG","SMUG1","KDM5B"), stack = T, flip = T, pt.size = 0, split.by = "sample") + NoLegend()

DotPlot(snscc_filtered, features = c("VEGFA","FLT1","NRP1","NRP2"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev)


gene_name <- c("S100A8")
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = gene_name, pt.size = 0.01, reduction = "umap.rna", colors_use = (brewer.pal(9,"Reds")), split.by = "sample")
#FeaturePlot(snscc_filtered, features = gene_name, pt.size = 0.01, reduction = "umap.wnn")
VlnPlot(snscc_filtered, features = gene_name, cols = color.snscc.filtered.CellType, pt.size = 0, split.by = "sample") + NoLegend()
#VlnPlot(snscc_filtered, features = gene_name, cols = color.snscc.filtered.CellType, pt.size = 0) + NoLegend()

#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("TP53","EGFR","CDKN2A","FGFR1","SOX2","NFKB1","PTGS2"), num_columns = 3, pt.size = 0.005, reduction = "umap.wnn", colors_use = rev(brewer.pal(9,"RdBu")))

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("CDH1","PDGFRB","CDH5","ITGAX","CD2","FCRL5"), num_columns = 3, pt.size = 0.005, reduction = "umap.wnn", colors_use = (brewer.pal(9,"Reds"))) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/snscc_filtered_res_12_featureplot_8*13.pdf", plot = figure, height = 8, width = 13, units = "in")

#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("EGFR","NDRG1","COL1A1","COL6A1","MSN","PNP","IFI30","TYROBP","IL32","SLC38A1"), num_columns = 2, pt.size = 0.005, reduction = "umap.wnn", colors_use = rev(brewer.pal(11,"RdYlBu")))-> figure
#figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_DEG_with_bulk_featureplot_8*18.pdf", plot = figure, height = 18, width = 8, units = "in")

# c("ADM","CALCRL")
# c("VEGFA","FLT1") FLT1=VEGFR1
# c("VEGFA","KDR") KDR=VEGFR2
# c("SPP1","ITGAV","ITGB1","ITGA5")
# c("CD44","RUNX1","NR2F2","DLL4","CALCRL","KDR")
DefaultAssay(snscc_filtered) <- "RNA"
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ADM","CALCRL","VEGFA","KDR"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = brewer.pal(n=9,"OrRd")) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_ADM_VEGFA_featureplot_9*8.pdf", plot = figure, width = 9, height = 8, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ADM","CALCRL","VEGFA","FLT1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#D7D8D9","#DDE1E2","#E8D8B4","#EFBB56","#F5A120","#DD7027","#C13424","#B01647","#9B217A"))
#"#E5E8EA"
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ADM","CALCRL","VEGFA","KDR"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#002153","#00366C","#0072B4","#79ABE2","#FFFFBF","#D07C42","#C7522B","#b21704","#611300")) -> figure
figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/snscc_filtered_res_12_ADM_VEGFA_featureplot_9*8.pdf", plot = figure, width = 9, height = 8, units = "in")



#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("SPP1","CD44","ITGAV","ITGB1","ITGA5"), num_columns = 5, pt.size = 0.005, reduction = "umap.wnn", colors_use = c("#D7D8D9","#DDE1E2","#E8D8B4","#EFBB56","#F5A120","#DD7027","#C13424","#B01647","#9B217A"))

#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NAMPT","INSR","ITGAV","ITGB1","ITGA5"), num_columns = 5, pt.size = 0.005, reduction = "umap.wnn", colors_use = c("#D7D8D9","#DDE1E2","#E8D8B4","#EFBB56","#F5A120","#DD7027","#C13424","#B01647","#9B217A"))

#FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("TGFB1","TGFBI","TGFBR1","TGFBR2","TGFBR3"), num_columns = 2, pt.size = 0.005, reduction = "umap.wnn", colors_use = c("#D7D8D9","#DDE1E2","#E8D8B4","#EFBB56","#F5A120","#DD7027","#C13424","#B01647","#9B217A"))

VlnPlot(snscc_filtered, features = c("ADM2"), pt.size = 0) + NoLegend() #+ UnRotate_X()

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NAMPT"), num_columns = 1, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"Reds")))


# DEGs (CellType)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")
snscc_filtered.all.genes <- FindAllMarkers(snscc_filtered, min.pct = 0.25, only.pos = T, logfc.threshold = 1)
snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05 & snscc_filtered.all.genes$avg_log2FC > 1, ]
snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered.da.genes.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.genes$gene, assays = "RNA")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(sample,CellType_sample) %>% unique()
anno_col
rownames(anno_col) <- anno_col$CellType_sample
anno_col
anno_colors <- list(CellType_sample = c("Epithelial cells(Tumor)" = "#D05353","Fibroblasts(Tumor)" = "#EA9014","Endothelial cells(Tumor)" = "#FCCF55","Myeloid cells(Tumor)" = "#76AB62","T cells(Tumor)" = "#4d648d","Plasma cells(Tumor)"="#9F8CB3","Epithelial cells(Normal)" = "#D05353","Fibroblasts(Normal)" = "#EA9014","Endothelial cells(Normal)" = "#FCCF55","Myeloid cells(Normal)" = "#76AB62","T cells(Normal)" = "#4d648d","Plasma cells(Normal)"="#9F8CB3"), sample = c("Normal"="#4A6274", "Tumor"="#E2725A"))
pheatmap::pheatmap(snscc_filtered.da.genes.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100)) -> figure
write.table(snscc_filtered.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_DEGs_log2FC_1.txt", sep = '\t')
#write.table(snscc_filtered.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DEGs_log2FC_1.txt", sep = '\t')
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_DEGs_Log2FC_1_heatmap_4*5.pdf", plot = figure, width = 4, height = 5, units = "in")

gene_counts <- table(snscc_filtered.all.da.genes$cluster)
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Gene number")

# DEGs (sample)
DefaultAssay(snscc_filtered) <- "RNA"
Idents(snscc_filtered) <- "sample"
snscc_filtered.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0, ident.1 = "Tumor", ident.2 = "HC")
write.table(snscc_filtered.all.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_sample_FindMarkers_DEGs_for_volcano_plot.txt", sep = '\t')

snscc_filtered.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Tumor", ident.2 = "HC")
snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05, ]
snscc_filtered.da.genes.average <- AverageExpression(snscc_filtered, features = rownames(snscc_filtered.all.da.genes), assays = "RNA")
pheatmap::pheatmap(snscc_filtered.da.genes.average[["RNA"]], scale = 'row', cluster_rows = T, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100))
write.table(snscc_filtered.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_sample_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

# DEGs (CellType_sample;FindAllMarkers)
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
#levels(snscc_filtered) <- c("Epithelial cells(HC)","Epithelial cells(Tumor)","Fibroblasts(HC)","Fibroblasts(Tumor)","Endothelial cells(HC)","Endothelial cells(Tumor)","Myeloid cells(HC)","Myeloid cells(Tumor)","T cells(HC)","T cells(Tumor)","Plasma cells(HC)","Plasma cells(Tumor)")
DefaultAssay(snscc_filtered) <- "RNA"
#snscc_filtered.all.genes <- FindAllMarkers(snscc_filtered, min.pct = 0.25, only.pos = T, logfc.threshold = 0.58)
#snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05, ]
#write.table(snscc_filtered.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DEGs_log2FC_058_HC_Tumor.txt", sep = '\t')

#snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05 & snscc_filtered.all.genes$avg_log2FC > 1, ]
#write.table(snscc_filtered.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DEGs_log2FC_1_HC_Tumor.txt", sep = '\t')

snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DEGs_log2FC_1_HC_Tumor.txt", sep = '\t')
snscc_filtered.da.genes.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.genes$gene, assays = "RNA")
pheatmap::pheatmap(snscc_filtered.da.genes.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_DEGs_Log2FC_058_heatmap_2*4.5.pdf", plot = figure, width = 2, height = 4.5, units = "in")
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_DEGs_Log2FC_1_heatmap_2*5.5.pdf", plot = figure, width = 2, height = 5.5, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_DEGs_Log2FC_1_heatmap_2*3.5.pdf", plot = figure, width = 2, height = 3.5, units = "in")

gene_counts <- table(snscc_filtered.all.da.genes$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType_sample", ylab = "Gene number")

# DEGs (CellType_sample;FindMarkers)
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Epithelial cells(Tumor)","Fibroblasts(HC)","Fibroblasts(Tumor)","Endothelial cells(HC)","Endothelial cells(Tumor)","Myeloid cells(HC)","Myeloid cells(Tumor)","T cells(HC)","T cells(Tumor)","Plasma cells(HC)","Plasma cells(Tumor)")
DefaultAssay(snscc_filtered) <- "RNA"

snscc_filtered_Epi.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Epithelial cells(Tumor)", ident.2 = "Epithelial cells(HC)")
snscc_filtered_Epi.all.da.genes <- snscc_filtered_Epi.all.genes[snscc_filtered_Epi.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Epi.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Epi_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_Fibro.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Fibroblasts(Tumor)", ident.2 = "Fibroblasts(HC)")
snscc_filtered_Fibro.all.da.genes <- snscc_filtered_Fibro.all.genes[snscc_filtered_Fibro.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Fibro.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_Endo.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Endothelial cells(Tumor)", ident.2 = "Endothelial cells(HC)")
snscc_filtered_Endo.all.da.genes <- snscc_filtered_Endo.all.genes[snscc_filtered_Endo.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Endo.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Endo_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_Mye.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Myeloid cells(Tumor)", ident.2 = "Myeloid cells(HC)")
snscc_filtered_Mye.all.da.genes <- snscc_filtered_Mye.all.genes[snscc_filtered_Mye.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Mye.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Mye_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_T.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "T cells(Tumor)", ident.2 = "T cells(HC)")
snscc_filtered_T.all.da.genes <- snscc_filtered_T.all.genes[snscc_filtered_T.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_T.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_T_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_PC.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "Plasma cells(Tumor)", ident.2 = "Plasma cells(HC)")
snscc_filtered_PC.all.da.genes <- snscc_filtered_PC.all.genes[snscc_filtered_PC.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_PC.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_PC_FindMarkers_DEGs_log2FC_058.txt", sep = '\t')

snscc_filtered_Epi.all.genes$cluster <- "Epi"
snscc_filtered_Fibro.all.genes$cluster <- "Fibro"
snscc_filtered_Endo.all.genes$cluster <- "Endo"
snscc_filtered_Mye.all.genes$cluster <- "Mye"
snscc_filtered_T.all.genes$cluster <- "T"
snscc_filtered_PC.all.genes$cluster <- "PC"

snscc_filtered_Epi.all.genes$gene <- rownames(snscc_filtered_Epi.all.genes)
snscc_filtered_Fibro.all.genes$gene <- rownames(snscc_filtered_Fibro.all.genes)
snscc_filtered_Endo.all.genes$gene <- rownames(snscc_filtered_Endo.all.genes)
snscc_filtered_Mye.all.genes$gene <- rownames(snscc_filtered_Mye.all.genes)
snscc_filtered_T.all.genes$gene <- rownames(snscc_filtered_T.all.genes)
snscc_filtered_PC.all.genes$gene <- rownames(snscc_filtered_PC.all.genes)


snscc_filtered_res_08_FindMarkers_DEGs_log2FC_058 <- rbind(snscc_filtered_Epi.all.genes,snscc_filtered_Fibro.all.genes,snscc_filtered_Endo.all.genes,snscc_filtered_Mye.all.genes,snscc_filtered_T.all.genes,snscc_filtered_PC.all.genes)
write.table(snscc_filtered_res_08_FindMarkers_DEGs_log2FC_058, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_058.txt", sep = "\t")


snscc_filtered_Epi.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0, ident.1 = "Epithelial cells(Tumor)", ident.2 = "Epithelial cells(HC)")
#snscc_filtered_Epi.all.da.genes <- snscc_filtered_Epi.all.genes[snscc_filtered_Epi.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Epi.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Epi_FindMarkers_DEGs_log2FC_0.txt", sep = '\t')

snscc_filtered_Fibro.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0, ident.1 = "Fibroblasts(Tumor)", ident.2 = "Fibroblasts(HC)")
#snscc_filtered_Fibro.all.da.genes <- snscc_filtered_Fibro.all.genes[snscc_filtered_Fibro.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Fibro.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_0.txt", sep = '\t')

snscc_filtered_Endo.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0, ident.1 = "Endothelial cells(Tumor)", ident.2 = "Endothelial cells(HC)")
#snscc_filtered_Endo.all.da.genes <- snscc_filtered_Endo.all.genes[snscc_filtered_Endo.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Endo.all.da.genes, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Endo_FindMarkers_DEGs_log2FC_0.txt", sep = '\t')

snscc_filtered_Epi.all.genes$cluster <- "Epi"
snscc_filtered_Fibro.all.genes$cluster <- "Fibro"
snscc_filtered_Endo.all.genes$cluster <- "Endo"

snscc_filtered_Epi.all.genes$gene <- rownames(snscc_filtered_Epi.all.genes)
snscc_filtered_Fibro.all.genes$gene <- rownames(snscc_filtered_Fibro.all.genes)
snscc_filtered_Endo.all.genes$gene <- rownames(snscc_filtered_Endo.all.genes)

snscc_filtered_res_08_FindMarkers_DEGs_log2FC_0 <- rbind(snscc_filtered_Epi.all.genes,snscc_filtered_Fibro.all.genes,snscc_filtered_Endo.all.genes)
write.table(snscc_filtered_res_08_FindMarkers_DEGs_log2FC_0, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_0_3CellType.txt", sep = "\t")

###############################################################






# DEGs (CellType_sample: heatmap)
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")

snscc_filtered_without_PC <- subset(snscc_filtered, CellType != "Plasma cells")
Idents(snscc_filtered_without_PC) <- "CellType"
levels(snscc_filtered_without_PC) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells")

snscc_filtered_Epi.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Epi_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered_Fibro.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered_Endo.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Endo_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered_Mye.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Mye_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered_T.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_T_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')
snscc_filtered_PC.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_PC_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_Epi.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Epi_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')
snscc_filtered_Fibro.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')
snscc_filtered_Endo.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Endo_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')
snscc_filtered_Mye.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_Mye_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')
snscc_filtered_T.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_T_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')
snscc_filtered_PC.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_15_CellType_sample_PC_FindMarkers_DEGs_log2FC_no_cutoff.txt", sep = '\t')


snscc_filtered_Epi.all.da.genes$cluster <- "Epithelial cells"
snscc_filtered_Epi.all.da.genes$Gene <- rownames(snscc_filtered_Epi.all.da.genes)
snscc_filtered_Epi.all.da.genes.order <- snscc_filtered_Epi.all.da.genes[order(snscc_filtered_Epi.all.da.genes$avg_log2FC), ]

snscc_filtered_Fibro.all.da.genes$cluster <- "Fibroblasts"
snscc_filtered_Fibro.all.da.genes$Gene <- rownames(snscc_filtered_Fibro.all.da.genes)
snscc_filtered_Fibro.all.da.genes.order <- snscc_filtered_Fibro.all.da.genes[order(snscc_filtered_Fibro.all.da.genes$avg_log2FC), ]

snscc_filtered_Endo.all.da.genes$cluster <- "Endothelial cells"
snscc_filtered_Endo.all.da.genes$Gene <- rownames(snscc_filtered_Endo.all.da.genes)
snscc_filtered_Endo.all.da.genes.order <- snscc_filtered_Endo.all.da.genes[order(snscc_filtered_Endo.all.da.genes$avg_log2FC), ]

snscc_filtered_Mye.all.da.genes$cluster <- "Myeloid cells"
snscc_filtered_Mye.all.da.genes$Gene <- rownames(snscc_filtered_Mye.all.da.genes)
snscc_filtered_Mye.all.da.genes.order <- snscc_filtered_Mye.all.da.genes[order(snscc_filtered_Mye.all.da.genes$avg_log2FC), ]

snscc_filtered_T.all.da.genes$cluster <- "T cells"
snscc_filtered_T.all.da.genes$Gene <- rownames(snscc_filtered_T.all.da.genes)
snscc_filtered_T.all.da.genes.order <- snscc_filtered_T.all.da.genes[order(snscc_filtered_T.all.da.genes$avg_log2FC), ]

snscc_filtered.CellType.FindMarkers <- rbind(snscc_filtered_Epi.all.da.genes.order,snscc_filtered_Fibro.all.da.genes.order,snscc_filtered_Endo.all.da.genes.order,snscc_filtered_Mye.all.da.genes.order,snscc_filtered_T.all.da.genes.order)
#write.table(snscc_filtered.CellType.FindMarkers, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_no_cutoff_rbind_without_PC_down_up.txt", sep = "\t")
#write.table(snscc_filtered.CellType.FindMarkers, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_058_rbind_without_PC_down_up.txt", sep = "\t")
#write.table(snscc_filtered.CellType.FindMarkers, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_rbind_without_PC_down_up.txt", sep = "\t")

snscc_filtered.CellType.FindMarkers.1 <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_rbind_without_PC_down_up.txt", sep = "\t")

snscc_filtered.CellType.FindMarkers <- snscc_filtered.CellType.FindMarkers[snscc_filtered.CellType.FindMarkers$Gene %in% snscc_filtered.CellType.FindMarkers.1$Gene, ]

#snscc_filtered_without_PC <- subset(snscc_filtered, CellType != "Plasma cells")

heatmap_data <- reshape(snscc_filtered.CellType.FindMarkers, idvar = "Gene", timevar = "cluster", direction = "wide")
heatmap_data <- heatmap_data[, c("Gene","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells","avg_log2FC.Myeloid cells","avg_log2FC.T cells")]
rownames(heatmap_data) <- heatmap_data$Gene

anno_col <- snscc_filtered_without_PC@meta.data %>% data.frame() %>% select(CellType) %>% unique()
anno_col
rownames(anno_col) <- c("avg_log2FC.Epithelial.cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial.cells","avg_log2FC.T.cells","avg_log2FC.Myeloid.cells")
anno_col
anno_colors <- list(CellType = c("Epithelial cells" = "#d42a34","Fibroblasts" = "#f78c37","Endothelial cells" = "#ffd827","Myeloid cells" = "#62ca50","T cells" = "#0677ba"))

heatmap_data <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_rbind_without_PC_down_up_heatmap_data_rearrange.txt", sep = "\t", header = T, row.names = 1)
#heatmap_data[is.na(heatmap_data)] <- 0

pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, na_col = "white", annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
#write.table(heatmap_data, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_058_rbind_without_PC_down_up_heatmap_data.txt", sep = "\t")
#write.table(heatmap_data, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_rbind_without_PC_down_up_heatmap_data.txt", sep = "\t")
#write.table(heatmap_data, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_no_cutoff_rbind_without_PC_down_up_heatmap_data.txt", sep = "\t")
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_058_down_up_heatmap_without_PC_3.5*6.pdf", plot = figure, width = 3.5, height = 6, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_down_up_heatmap_without_PC_3.5*6.pdf", plot = figure, width = 3.5, height = 6, units = "in")

gene_counts <- table(snscc_filtered.all.da.genes$Concat)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Gene number")

# DARs (CellType_sample: heatmap)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")

DefaultAssay(snscc_filtered_without_PC) <- "peaks"
Idents(snscc_filtered_without_PC) <- "CellType"
levels(snscc_filtered_without_PC) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells")

snscc_filtered_Epi.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_058_minpct_001.txt", sep = '\t')
snscc_filtered_Fibro.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_058_minpct_001.txt", sep = '\t')
snscc_filtered_Endo.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_058_minpct_001.txt", sep = '\t')
snscc_filtered_Mye.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_058_minpct_001.txt", sep = '\t')
snscc_filtered_T.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_058_minpct_001.txt", sep = '\t')

snscc_filtered_Epi.all.da.peaks$cluster <- "Epithelial cells"
snscc_filtered_Epi.all.da.peaks$Gene <- rownames(snscc_filtered_Epi.all.da.peaks)
snscc_filtered_Epi.all.da.peaks.order <- snscc_filtered_Epi.all.da.peaks[order(-snscc_filtered_Epi.all.da.peaks$avg_log2FC), ]

snscc_filtered_Fibro.all.da.peaks$cluster <- "Fibroblasts"
snscc_filtered_Fibro.all.da.peaks$Gene <- rownames(snscc_filtered_Fibro.all.da.peaks)
snscc_filtered_Fibro.all.da.peaks.order <- snscc_filtered_Fibro.all.da.peaks[order(-snscc_filtered_Fibro.all.da.peaks$avg_log2FC), ]

snscc_filtered_Endo.all.da.peaks$cluster <- "Endothelial cells"
snscc_filtered_Endo.all.da.peaks$Gene <- rownames(snscc_filtered_Endo.all.da.peaks)
snscc_filtered_Endo.all.da.peaks.order <- snscc_filtered_Endo.all.da.peaks[order(-snscc_filtered_Endo.all.da.peaks$avg_log2FC), ]

snscc_filtered_Mye.all.da.peaks$cluster <- "Myeloid cells"
snscc_filtered_Mye.all.da.peaks$Gene <- rownames(snscc_filtered_Mye.all.da.peaks)
snscc_filtered_Mye.all.da.peaks.order <- snscc_filtered_Mye.all.da.peaks[order(-snscc_filtered_Mye.all.da.peaks$avg_log2FC), ]

snscc_filtered_T.all.da.peaks$cluster <- "T cells"
snscc_filtered_T.all.da.peaks$Gene <- rownames(snscc_filtered_T.all.da.peaks)
snscc_filtered_T.all.da.peaks.order <- snscc_filtered_T.all.da.peaks[order(-snscc_filtered_T.all.da.peaks$avg_log2FC), ]

snscc_filtered.CellType.FindMarkers.peaks <- rbind(snscc_filtered_Epi.all.da.peaks.order,snscc_filtered_Fibro.all.da.peaks.order,snscc_filtered_Endo.all.da.peaks.order,snscc_filtered_Mye.all.da.peaks.order,snscc_filtered_T.all.da.peaks.order) # plasma cell DAR number = 0
write.table(snscc_filtered.CellType.FindMarkers.peaks, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DARs_log2FC_058_rbind_minpct_001_up_down_without_PC.txt", sep = "\t")
#write.table(snscc_filtered.CellType.FindMarkers.peaks, "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DARs_log2FC_1_rbind.txt", sep = "\t")

#snscc_filtered.CellType.FindMarkers.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1_rbind.txt", sep = "\t")

snscc_filtered.CellType.FindMarkers.peaks.up <- snscc_filtered.CellType.FindMarkers.peaks[snscc_filtered.CellType.FindMarkers.peaks$avg_log2FC > 0, ]

snscc_filtered.CellType.FindMarkers.peaks.down <- snscc_filtered.CellType.FindMarkers.peaks[snscc_filtered.CellType.FindMarkers.peaks$avg_log2FC < 0, ]

heatmap_data <- reshape(snscc_filtered.CellType.FindMarkers.peaks, idvar = "Gene", timevar = "cluster", direction = "wide")
heatmap_data <- heatmap_data[, c("Gene","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells","avg_log2FC.Myeloid cells","avg_log2FC.T cells")]
rownames(heatmap_data) <- heatmap_data$Gene
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
anno_col
rownames(anno_col) <- c("avg_log2FC.Fibroblasts","avg_log2FC.Epithelial cells","avg_log2FC.Endothelial cells","avg_log2FC.T cells","avg_log2FC.Myeloid cells")
anno_col
anno_colors <- list(CellType = c("Epithelial cells" = "#d42a34","Fibroblasts" = "#f78c37","Endothelial cells" = "#ffd827","Myeloid cells" = "#62ca50","T cells" = "#0677ba"))

pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, na_col = "white", annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_08_CellType_sample_FindMarkers_DARs_Log2FC_1_heatmap_3*5.pdf", plot = figure, width = 3.5, height = 6, units = "in")

gene_counts <- table(snscc_filtered.CellType.FindMarkers.peaks.up$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Gene number")

###
### Correlation heatmap
Idents(snscc_filtered) <- "CellType_sample"
snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DEGs_log2FC_1_HC_Tumor.txt", sep = '\t')
snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_DEGs_log2FC_1.txt", sep = '\t')
DefaultAssay(snscc_filtered) <- "RNA"
avg_snscc_filtered.all.da.genes <- AverageExpression(snscc_filtered, assays = "RNA", slot = "data", features = intersect(snscc_filtered.all.da.genes$gene, rownames(snscc_filtered)))

DefaultAssay(snscc_filtered) <- "RNA"
avg_snscc_filtered.all.da.genes <- AverageExpression(snscc_filtered, assays = "RNA", slot = "data", features = VariableFeatures(snscc_filtered))

#VariableFeatures is original code instead of rownames because peak number in VariableFeatures and rownames is same in peaks assay
cor_mat_snscc_filtered.all.genes <- cor(avg_snscc_filtered.all.da.genes$RNA, method = "spearman")
#anno_col <- Epi@meta.data %>% data.frame() %>% select(seurat_clusters) %>% unique()
#rownames(anno_col) <- anno_col$seurat_clusters
#anno_colors <- list(seurat_clusters = c("0" = "#C79CA6","1" = "#7D93B5","2" = "#ECC371","3" = "#C4BB33","4" = "#D05353","5"="#8EB18F","6"="#1b3461","7"="#bebada","8"="#B3906C","9"="#C7B8A4"))
pheatmap::pheatmap(cor_mat_snscc_filtered.all.genes, cluster_rows = T, cluster_cols = T, show_rownames = T, treeheight_row = 30, treeheight_col = 30, border_color = "black", display_numbers = T, number_color = "black", color = colorRampPalette(rev(brewer.pal(n=10,name="RdYlBu")))(100))

# DARs (CellType)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType"
#snscc_filtered.all.peaks <- FindAllMarkers(snscc_filtered, min.pct = 0.02, only.pos = T, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks')
#snscc_filtered.all.da.peaks <- snscc_filtered.all.peaks[snscc_filtered.all.peaks$p_val_adj < 0.05 & snscc_filtered.all.peaks$avg_log2FC > 0.58, ]
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DARs_Log2FC_058.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
rownames(anno_col) <- anno_col$CellType
anno_colors <- list(CellType = c("Epithelial cells" = "#D05353","Fibroblasts" = "#EA9014","Endothelial cells" = "#FCCF55","Myeloid cells" = "#76AB62","T cells" = "#4d648d","Plasma cells" = "#9F8CB3"))
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_DARs_log2FC_058_heatmap_3.5*6.pdf", plot = figure, height = 6, width = 3.5, units = "in")
#write.table(snscc_filtered.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DARs_Log2FC_058.txt", sep = '\t')
da_peaks <- snscc_filtered.all.da.peaks$gene
peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
head(peaks_annotation)
write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DARs_Log2FC_058_annotation.txt", sep = '\t')

gene_counts <- table(snscc_filtered.all.da.peaks$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Peak number")

# reverse heatmap
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Plasma cells(Tumor)","Plasma cells(Normal)","T cells(Tumor)","T cells(Normal)","Myeloid cells(Tumor)","Myeloid cells(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Epithelial cells(Tumor)","Epithelial cells(Normal)")
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_DARs_Log2FC_058.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(sample,CellType_sample) %>% unique()
anno_col
rownames(anno_col) <- anno_col$CellType_sample
anno_col
anno_colors <- list(CellType_sample = c("Epithelial cells(Tumor)" = "#D05353","Fibroblasts(Tumor)" = "#EA9014","Endothelial cells(Tumor)" = "#FCCF55","Myeloid cells(Tumor)" = "#76AB62","T cells(Tumor)" = "#4d648d","Plasma cells(Tumor)"="#9F8CB3","Epithelial cells(Normal)" = "#D05353","Fibroblasts(Normal)" = "#EA9014","Endothelial cells(Normal)" = "#FCCF55","Myeloid cells(Normal)" = "#76AB62","T cells(Normal)" = "#4d648d","Plasma cells(Normal)"="#9F8CB3"), sample = c("Normal"="#4A6274", "Tumor"="#E2725A"))
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
#levels(snscc_filtered) <- c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)","Myeloid cells(Tumor)","Myeloid cells(Normal)","T cells(Tumor)","T cells(Normal)","Plasma cells(Tumor)","Plasma cells(Normal)")

# DARs (CellType_sample)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

# TSS / Distal
## split peaks between promoter regions and distal elements in each CellType
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
tss.positions <- GetTSSPositions(ranges = annotation)
tss.positions <- Extend(x=tss.positions, upstream = 1000, downstream = 100, from.midpoint = T) # obtain region around TSS

snscc_filtered_Epi.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Epithelial cells(Tumor)", ident.2 = "Epithelial cells(HC)")
snscc_filtered_Epi.all.da.peaks <- snscc_filtered_Epi.all.peaks[snscc_filtered_Epi.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Epi.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_Epi.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_Epi_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_epi <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Epi.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_epi, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Fibro.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Fibroblasts(Tumor)", ident.2 = "Fibroblasts(HC)")
snscc_filtered_Fibro.all.da.peaks <- snscc_filtered_Fibro.all.peaks[snscc_filtered_Fibro.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Fibro.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_Fibro.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_Fibro_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_fibro <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Fibro.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_fibro, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Endo.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Endothelial cells(Tumor)", ident.2 = "Endothelial cells(HC)")
snscc_filtered_Endo.all.da.peaks <- snscc_filtered_Endo.all.peaks[snscc_filtered_Endo.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Endo.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_Endo.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_Endo_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_endo <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Endo.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_endo, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Mye.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Myeloid cells(Tumor)", ident.2 = "Myeloid cells(HC)")
snscc_filtered_Mye.all.da.peaks <- snscc_filtered_Mye.all.peaks[snscc_filtered_Mye.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Mye.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_Mye.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_Mye_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_mye <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Mye.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_mye, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_T.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "T cells(Tumor)", ident.2 = "T cells(HC)")
snscc_filtered_T.all.da.peaks <- snscc_filtered_T.all.peaks[snscc_filtered_T.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_T.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_T.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_T_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_t <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_T.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_t, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_PC.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Plasma cells(Tumor)", ident.2 = "Plasma cells(HC)")
snscc_filtered_PC.all.da.peaks <- snscc_filtered_PC.all.peaks[snscc_filtered_PC.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_PC.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_PC_DARs_log2FC_058_minpct_002.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered_PC.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_PC_DARs_log2FC_058_annotation.txt", sep = '\t')
closest_tss_snscc_filtered_pc <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_PC.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_pc, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_PC_DARs_log2FC_058_minpct_002_closest_tss.txt", sep = "\t")

# DARs (CellType_sample: heatmap/log2FC)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")
snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_FindMarkers_DARs_log2FC_1_minpct_002_for_log_heatmap.txt", sep = "\t", header = T)
Epi.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "Epithelial cells", ]
Epi.DARs <- Epi.DARs[order(Epi.DARs[, "avg_log2FC"], decreasing = T), ]
Fibro.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "Fibroblasts", ]
Fibro.DARs <- Fibro.DARs[order(Fibro.DARs[, "avg_log2FC"], decreasing = T), ]
Endo.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "Endothelial cells", ]
Endo.DARs <- Endo.DARs[order(Endo.DARs[, "avg_log2FC"], decreasing = T), ]
Mye.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "Myeloid cells", ]
Mye.DARs <- Mye.DARs[order(Mye.DARs[, "avg_log2FC"], decreasing = F), ]
T.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "T cells", ]
T.DARs <- T.DARs[order(T.DARs[, "avg_log2FC"], decreasing = F), ]
PC.DARs <- snscc_filtered.all.da.genes[snscc_filtered.all.da.genes$cluster == "Plasma cells", ]
PC.DARs <- PC.DARs[order(PC.DARs[, "avg_log2FC"], decreasing = F), ]
PC.DARs <- rbind(PC.DARs, numeric(ncol(PC.DARs)))
rownames(PC.DARs)[nrow(PC.DARs)] <- "chrY-7273556-7274399"
rownames(PC.DARs)
colnames(PC.DARs) <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","Region")
PC.DARs$cluster <- "Plasma cells"
PC.DARs$Region <- "chrY-7273556-7274399"

#snscc_filtered.all.da.genes.re <- rbind(Epi.DARs,Fibro.DARs,Endo.DARs,Mye.DARs,T.DARs,PC.DARs)
snscc_filtered.all.da.genes.re <- rbind(Epi.DARs,Fibro.DARs,Endo.DARs)

heatmap_data <- reshape(snscc_filtered.all.da.genes.re, idvar = "Region", timevar = "cluster", direction = "wide")
#heatmap_data <- heatmap_data[, c("Region","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells","avg_log2FC.Myeloid cells","avg_log2FC.T cells","avg_log2FC.Plasma cells")]
heatmap_data <- heatmap_data[, c("Region","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells")]
rownames(heatmap_data) <- heatmap_data$Region
#anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, fontsize_row = 5, na_col = "white", color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_FindMarkers_DARs_Log2FC_1_minpct_002_3CellType_heatmap_4*6.pdf", plot = figure, width = 4, height = 6, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_FindMarkers_DARs_Log2FC_1_minpct_002_All_CellType_heatmap_8*8.pdf", plot = figure, width = 8, height = 8, units = "in")

snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_FindMarkers_DARs_log2FC_1_minpct_002.txt", sep = "\t", header = T)
gene_counts <- table(snscc_filtered.all.da.genes$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Gene number")

heatmap_data <- reshape(T.DARs, idvar = "Region", timevar = "cluster", direction = "wide")
heatmap_data <- heatmap_data[, c("Region","avg_log2FC.T cells")]
pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, na_col = "white", color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure

# DEGs (CellType_sample: heatmap/log2FC)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")
snscc_filtered.all.da.genes <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_FindMarkers_DARs_log2FC_1_heatmap.txt", sep = "\t", header = T)

heatmap_data <- reshape(snscc_filtered.all.da.genes, idvar = "Region", timevar = "cluster", direction = "wide")
heatmap_data <- heatmap_data[, c("Region","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells","avg_log2FC.Myeloid cells","avg_log2FC.T cells","avg_log2FC.Plasma cells")]
rownames(heatmap_data) <- heatmap_data$Region
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
anno_col
rownames(anno_col) <- c("avg_log2FC.Fibroblasts","avg_log2FC.Epithelial cells","avg_log2FC.Plasma cells","avg_log2FC.Endothelial cells","avg_log2FC.T cells","avg_log2FC.Myeloid cells")
anno_col
anno_colors <- list(CellType = c("Epithelial cells" = "#D05353","Fibroblasts" = "#EA9014","Endothelial cells" = "#FCCF55","Myeloid cells" = "#76AB62","T cells" = "#4d648d","Plasma cells"="#9F8CB3"))
pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, na_col = "white", annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_DEGs_LogFC_058_heatmap_3*5.pdf", plot = figure, width = 3, height = 5, units = "in")

gene_counts <- table(snscc_filtered.all.da.genes$cluster)
gene_counts
barplot(gene_counts, main = "Gene Counts per CellType", xlab = "CellType", ylab = "Gene number")

# DARs_FindAllMarkers_rearrange_HC_Tumor
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
snscc_filtered.all.peaks2 <- FindAllMarkers(snscc_filtered, min.pct = 0.01, only.pos = T, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks')
snscc_filtered.all.da.peaks2 <- snscc_filtered.all.peaks2[snscc_filtered.all.peaks2$p_val_adj < 0.05 & snscc_filtered.all.peaks2$avg_log2FC > 0.58, ]
write.table(snscc_filtered.all.da.peaks2, file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DARs_Log2FC_058_minpct_001_HC_Tumor.txt", sep = '\t')
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Documents/snscc_filtered_res_08_CellType_sample_DARs_Log2FC_058_minpct_002_HC_Tumor.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="BrBG")))(100)) -> figure
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(ArchRPalettes$blueYellow)(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_rearrange_HC_Tumor_heatmap_2*6.png", plot = figure, height = 6, width = 2, units = "in")
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_rearrange_HC_Tumor_heatmap_3*7.pdf", plot = figure, height = 7, width = 3, units = "in")
#write.table(snscc_filtered.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_rearrange_HC_Tumor.txt", sep = '\t')
# peak annotation
da_peaks <- snscc_filtered.all.da.peaks$gene
peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_rearrange_HC_Tumor_annotation.txt", sep = '\t')

# FindAllMarkers
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")
#snscc_filtered.all.peaks <- FindAllMarkers(snscc_filtered, min.pct = 0.02, only.pos = T, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks')
#snscc_filtered.all.da.peaks <- snscc_filtered.all.peaks[snscc_filtered.all.peaks$p_val_adj < 0.05 & snscc_filtered.all.peaks$avg_log2FC > 0.58, ]
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(sample,CellType_sample) %>% unique()
anno_col
rownames(anno_col) <- anno_col$CellType_sample
anno_col
anno_colors <- list(CellType_sample = c("Epithelial cells(Tumor)" = "#D05353","Fibroblasts(Tumor)" = "#EA9014","Endothelial cells(Tumor)" = "#FCCF55","Myeloid cells(Tumor)" = "#76AB62","T cells(Tumor)" = "#4d648d","Plasma cells(Tumor)"="#9F8CB3","Epithelial cells(Normal)" = "#D05353","Fibroblasts(Normal)" = "#EA9014","Endothelial cells(Normal)" = "#FCCF55","Myeloid cells(Normal)" = "#76AB62","T cells(Normal)" = "#4d648d","Plasma cells(Normal)"="#9F8CB3"), sample = c("Normal"="#4A6274", "Tumor"="#E2725A"))
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
#pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
write.table(snscc_filtered.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse.txt", sep = '\t')
# peak annotation
da_peaks <- snscc_filtered.all.da.peaks$gene
peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse_annotation.txt", sep = '\t')
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_heatmap_4.5*6.png", plot = figure, height = 6, width = 4.5, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_heatmap_6*10.pdf", plot = figure, height = 10, width = 6, units = "in")
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_TSS_logFC_058_minpct_002_reverse_heatmap_4.5*3.pdf", plot = figure, height = 3, width = 4.5, units = "in")

# reverse heatmap
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Plasma cells(Tumor)","Plasma cells(Normal)","T cells(Tumor)","T cells(Normal)","Myeloid cells(Tumor)","Myeloid cells(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Epithelial cells(Tumor)","Epithelial cells(Normal)")
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_sample_DARs_log2FC_058_heatmap_reverse_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(sample,CellType_sample) %>% unique()
anno_col
rownames(anno_col) <- anno_col$CellType_sample
anno_col
anno_colors <- list(CellType_sample = c("Epithelial cells(Tumor)" = "#D05353","Fibroblasts(Tumor)" = "#EA9014","Endothelial cells(Tumor)" = "#FCCF55","Myeloid cells(Tumor)" = "#76AB62","T cells(Tumor)" = "#4d648d","Plasma cells(Tumor)"="#9F8CB3","Epithelial cells(Normal)" = "#D05353","Fibroblasts(Normal)" = "#EA9014","Endothelial cells(Normal)" = "#FCCF55","Myeloid cells(Normal)" = "#76AB62","T cells(Normal)" = "#4d648d","Plasma cells(Normal)"="#9F8CB3"), sample = c("Normal"="#4A6274", "Tumor"="#E2725A"))
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_heatmap_reverse_4.5*6.png", plot = figure, height = 6, width = 4.5, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_DEPs_logFC_058_minpct_002_heatmap_reverse_6*10.pdf", plot = figure, height = 10, width = 6, units = "in")
gene_counts <- table(snscc_filtered.all.da.peaks$cluster)
barplot(gene_counts, main = "Gene Counts per Cluster", xlab = "CellType", ylab = "Peak number")
levels(snscc_filtered) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")

# TSS / Distal
## split peaks between promoter regions and distal elements in each CellType
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
library(stringr)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
ucsc.levels
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"
genome(annotation)
seqlevelsStyle(annotation)

DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
tss.positions <- GetTSSPositions(ranges = annotation)
tss.positions <- Extend(x=tss.positions, upstream = 1000, downstream = 100, from.midpoint = T) # obtain region around TSS

DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse.txt", sep = '\t')
snscc_filtered_DEPs <- snscc_filtered.all.da.peaks$gene
snscc_filtered_DEPs
closest_tss_snscc_filtered <- ClosestFeature(snscc_filtered, snscc_filtered_DEPs, annotation = tss.positions)
write.table(closest_tss_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse_closest_tss.txt", sep = "\t")

DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002.txt", sep = '\t')
snscc_filtered_DEPs <- snscc_filtered.all.da.peaks$gene
snscc_filtered_DEPs
closest_tss_snscc_filtered <- ClosestFeature(snscc_filtered, snscc_filtered_DEPs, annotation = tss.positions)
write.table(closest_tss_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_CellType_sample_DARs_minpct_002_logFC_058_closest_tss.txt", sep = "\t")
library(dplyr)
tss_peaks_snscc_filtered <- closest_tss_snscc_filtered %>% filter(distance == 0) # peaks falling within promoter region (between -1000 and +100 from TSS)
write.table(tss_peaks_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_CellType_sample_DARs_minpct_002_logFC_058_tss_peaks.txt", sep = "\t")
distal_peaks_snscc_filtered <- closest_tss_snscc_filtered %>% filter(distance > 0) # peaks falling outise promoter region - assigned as distal
write.table(distal_peaks_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_CellType_sample_DARs_minpcts_002_logFC_058_distal_peaks.txt", sep = "\t")

# DARs (sample)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "sample"
#snscc_filtered.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 0.58, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Tumor", ident.2 = "Normal")
#snscc_filtered.all.da.peaks <- snscc_filtered.all.peaks[snscc_filtered.all.peaks$p_val_adj < 0.05, ]
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_sample_DARs_Log2FC_058.txt", sep = '\t')
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene, assays = "peaks")
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
#write.table(snscc_filtered.all.da.peaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_sample_DARs_Log2FC_058.txt", sep = '\t')
# peak annotation
#da_peaks <- rownames(snscc_filtered.all.da.peaks)
#peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_sample_DARs_Log2FC_058_annotation.txt", sep = '\t')

# TSS / Distal
## split peaks between promoter regions and distal elements in each CellType
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
library(stringr)
ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
ucsc.levels
seqlevels(annotation) <- ucsc.levels
genome(annotation) <- "hg38"
genome(annotation)
seqlevelsStyle(annotation)

DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "sample"
tss.positions <- GetTSSPositions(ranges = annotation)
tss.positions <- Extend(x=tss.positions, upstream = 1000, downstream = 100, from.midpoint = T) # obtain region around TSS

DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered.all.da.peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_sample_DARs_Log2FC_058.txt", sep = '\t')
snscc_filtered_DEPs <- rownames(snscc_filtered.all.da.peaks)
snscc_filtered_DEPs
closest_tss_snscc_filtered <- ClosestFeature(snscc_filtered, snscc_filtered_DEPs, annotation = tss.positions)
write.table(closest_tss_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_sample_DARs_minpct_002_logFC_058_closest_tss.txt", sep = "\t")
library(dplyr)
tss_peaks_snscc_filtered <- closest_tss_snscc_filtered %>% filter(distance == 0) # peaks falling within promoter region (between -1000 and +100 from TSS)
write.table(tss_peaks_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_sample_DARs_minpct_002_logFC_058_tss_peaks.txt", sep = "\t")
distal_peaks_snscc_filtered <- closest_tss_snscc_filtered %>% filter(distance > 0) # peaks falling outise promoter region - assigned as distal
write.table(distal_peaks_snscc_filtered, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_sample_DARs_minpcts_002_logFC_058_distal_peaks.txt", sep = "\t")

# DEPs overlap with DEGs
library(pheatmap)
library(RColorBrewer)
library(dplyr)
DefaultAssay(snscc_filtered) <- "peaks"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells")
snscc_filtered.all.peaks <- read.table(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/snscc_filtered_CellType_DEPs_overlap_DEGs_minpct_002.txt", sep = '\t', header = T)
snscc_filtered.all.da.peaks <- snscc_filtered.all.peaks[snscc_filtered.all.peaks$p_val_adj < 0.05 & snscc_filtered.all.peaks$avg_log2FC > 0.58, ]
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$region, assays = "peaks")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
rownames(anno_col) <- anno_col$CellType
anno_colors <- list(CellType = c("Epithelial cells" = "#9F8CB3","Fibroblasts" = "#4d648d","Endothelial cells" = "#EA9014","Myeloid cells" = "#D05353","T cells" = "#76AB62"))
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, annotation_col = anno_col, annotation_colors = anno_colors, labels_row = snscc_filtered.all.da.peaks$region, color = colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_overlap_DEGs_logFC_058_minpct_002_peaks_2.6*4.pdf", plot = figure, height = 4, width = 2.6, units = "in")
DefaultAssay(snscc_filtered) <- "RNA"
snscc_filtered.da.peaks.average <- AverageExpression(snscc_filtered, features = snscc_filtered.all.da.peaks$gene_name, assays = "RNA")
pheatmap::pheatmap(snscc_filtered.da.peaks.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, annotation_col = anno_col, annotation_colors = anno_colors, labels_row = snscc_filtered.all.da.peaks$gene_name, color = colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_overlap_DEGs_logFC_058_minpct_002_rna_2.6*4.pdf", plot = figure, height = 4, width = 2.6, units = "in")

# peak annotation
da_peaks <- snscc_filtered.all.da.peaks$gene
peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
head(peaks_annotation)
write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/Epi/snscc_filtered_CellType_sample_DEPs_Log2FC_058_minpct_002_annotation.txt", sep = '\t')

# DEPs overlap with DEGs (multiome data)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058 <- read.table(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.txt", header = T, sep = '\t')
snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.average <- AverageExpression(snscc_filtered, features = snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058$gene, assays = "peaks")
anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(CellType) %>% unique()
rownames(anno_col) <- anno_col$CellType
anno_colors <- list(CellType = c("Epithelial cells" = "#9F8CB3","Fibroblasts" = "#4d648d","Endothelial cells" = "#EA9014","Myeloid cells" = "#D05353","T cells" = "#76AB62"))
pheatmap::pheatmap(snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058_heatmap_3*5.pdf", plot = figure, height = 5, width = 3, units = "in")
pheatmap::pheatmap(snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", fontsize_row = 5, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=7,name="RdYlBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058_heatmap_7*8.pdf", plot = figure, height = 8, width = 7, units = "in")

DefaultAssay(snscc_filtered) <- "RNA"
snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058 <- read.table(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.txt", header = T, sep = '\t')
snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.average <- AverageExpression(snscc_filtered, features = snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058$gene_name, assays = "RNA")
pheatmap::pheatmap(snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, border_color = "black", fontsize_row = 5, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)) -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058_RNA_assay_heatmap_3*5.pdf", plot = figure, height = 5, width = 3, units = "in")


## Upset Plot 
# https://zenn.dev/rchiji/articles/0cdf3fec190c31
DefaultAssay(snscc_filtered) <- "RNA"
snscc_filtered.all.genes <- read.table(file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DEGs_log2FC_1.txt", sep = '\t')#
#snscc_filtered.all.genes <- read.table(file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DARs_Log2FC_058_minpct_002_reverse.txt", sep = '\t')
DEG_list <- split(snscc_filtered.all.genes$gene, snscc_filtered.all.genes$cluster)
list_to_matrix(DEG_list)
library(ComplexHeatmap)
DEG_mt <- ComplexHeatmap::make_comb_mat(DEG_list)
DEG_mt
ComplexHeatmap::UpSet(DEG_mt)
#pdf(file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/snscc_filtered_CellType_DEGs_FindMarkers_UpsetPlot_3*16.pdf", height = 4, width = 16)
ComplexHeatmap::UpSet(DEG_mt, set_order = c("Epithelial cells_T","Fibroblasts_T","Endothelial cells_T","Myeloid cells_T","T cells_T","Plasma cells_T","Epithelial cells_N","Fibroblasts_N","Endothelial cells_N","Myeloid cells_N","T cells_N","Plasma cells_N"), row_names_gp=gpar(fontsize=10), right_annotation = upset_right_annotation(DEG_mt, gp=gpar(fill=c("#D05353","#EA9014","#FCCF55","#76AB62","#4d648d","#9F8CB3","#D05353","#EA9014","#FCCF55","#76AB62","#4d648d","#9F8CB3")), add_numbers=T, width = unit(2,"cm"), annotation_name_gp=gpar(fontsize=10)), top_annotation = upset_top_annotation(DEG_mt, add_numbers=T, annotation_name_gp=gpar(fontsize=10)))
dev.off()
# ggsave is not working


DefaultAssay(snscc_filtered) <- "RNA"
snscc_filtered.all.genes <- read.table(file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_CellType_sample_DEGs_log2FC_1_reverse_with_bulk_up_DEGs.txt", header = T, sep = '\t')
snscc_filtered.all.da.genes <- snscc_filtered.all.genes[snscc_filtered.all.genes$p_val_adj < 0.05 & snscc_filtered.all.genes$avg_log2FC > 1, ]
DEG_list <- split(snscc_filtered.all.genes$gene, snscc_filtered.all.genes$cluster)
list_to_matrix(DEG_list)
library(ComplexHeatmap)
DEG_mt <- ComplexHeatmap::make_comb_mat(DEG_list)
DEG_mt
ComplexHeatmap::UpSet(DEG_mt, add_numbers=T)
ComplexHeatmap:: UpSet(DEG_mt, set_order = c("Fibroblasts(Normal)","Endothelial cells(Normal)","Myeloid cells(Normal)","T cells(Normal)","Plasma cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)"), row_names_gp=gpar(fontsize=10), right_annotation = upset_right_annotation(DEG_mt, gp=gpar(fill=c("#D05353","#EA9014","#FCCF55","#76AB62","#4d648d","#9F8CB3","#D05353","#EA9014")), add_numbers=T, width = unit(3,"cm"), annotation_name_gp=gpar(fontsize=10)), top_annotation = upset_top_annotation(DEG_mt, add_numbers=T, annotation_name_gp=gpar(fontsize=10)))
pdf(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEGs_with_bulk_UpsetPlot_3*4.pdf", height = 4, width = 3)
ComplexHeatmap:: UpSet(t(DEG_mt), set_order = c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells"), row_names_gp=gpar(fontsize=10), top_annotation = upset_top_annotation(t(DEG_mt), gp=gpar(fill=c("#EA9014","#9F8CB3","#4d648d","#D05353","#76AB62")), add_numbers=F, width = unit(3,"cm"), annotation_name_gp=gpar(fontsize=10)), right_annotation = upset_right_annotation(t(DEG_mt), add_numbers=T, width = unit(3,"cm"), annotation_name_gp=gpar(fontsize=10)))
# ggsave is not working
dev.off()

DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered.all.peaks <- read.table(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/snscc_filtered_CellType_DEPs_LogFC_02_minpct_001.txt", sep = '\t')
snscc_filtered.all.da.peaks <- snscc_filtered.all.peaks[snscc_filtered.all.peaks$p_val_adj < 0.05 & snscc_filtered.all.peaks$avg_log2FC > 0.58, ]
DEP_list <- split(snscc_filtered.all.da.peaks$gene, snscc_filtered.all.da.peaks$cluster)
list_to_matrix(DEP_list)
library(ComplexHeatmap)
DEP_mt <- ComplexHeatmap::make_comb_mat(DEP_list)
DEP_mt
ComplexHeatmap::UpSet(DEP_mt)
ComplexHeatmap:: UpSet(DEP_mt, set_order = c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells"), row_names_gp=gpar(fontsize=10), right_annotation = upset_right_annotation(DEP_mt, gp=gpar(fill=c("#EA9014","#9F8CB3","#4d648d","#D05353","#76AB62")), add_numbers=T, width = unit(3.5,"cm"), annotation_name_gp=gpar(fontsize=10)), top_annotation = upset_top_annotation(DEP_mt, add_numbers=T, annotation_name_gp=gpar(fontsize=10)))

DefaultAssay(snscc_filtered) <- "peaks"
open.peaks <- read.table(file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open_peaks_CellType.txt", sep = '\t', header = T)
open_peak_list <- split(open.peaks$gene, open.peaks$cluster)
list_to_matrix(open_peak_list)
library(ComplexHeatmap)
open_peak_mt <- ComplexHeatmap::make_comb_mat(open_peak_list)
open_peak_mt
ComplexHeatmap::UpSet(open_peak_mt)
ComplexHeatmap:: UpSet(open_peak_mt, set_order = c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells"), row_names_gp=gpar(fontsize=10), right_annotation = upset_right_annotation(open_peak_mt, gp=gpar(fill=c("#EA9014","#9F8CB3","#4d648d","#D05353","#76AB62")), add_numbers=T, width = unit(2,"cm"), annotation_name_gp=gpar(fontsize=10)), top_annotation = upset_top_annotation(open_peak_mt, add_numbers=T, annotation_name_gp=gpar(fontsize=10)))


## Find accessible peaks in a set of cells
DefaultAssay(snscc_filtered) <- "peaks"
open.peaks.epi <- AccessiblePeaks(snscc_filtered, idents = "Epithelial cells")
open.peaks.fibro <- AccessiblePeaks(snscc_filtered, idents = "Fibroblasts")
open.peaks.endo <- AccessiblePeaks(snscc_filtered, idents = "Endothelial cells")
open.peaks.mye <- AccessiblePeaks(snscc_filtered, idents = "Myeloid cells")
open.peaks.t <- AccessiblePeaks(snscc_filtered, idents = "T cells")

write.table(open.peaks.epi, file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open.peaks.epi.txt", sep = "\t")
write.table(open.peaks.fibro, file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open.peaks.fibro.txt", sep = "\t")
write.table(open.peaks.endo, file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open.peaks.endo.txt", sep = "\t")
write.table(open.peaks.mye, file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open.peaks.mye.txt", sep = "\t")
write.table(open.peaks.t, file = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/update/open.peaks.t.txt", sep = "\t")


## Linking peaks to genes
DefaultAssay(snscc_filtered) <- "peaks"
#Idents(snscc_filtered) <- "CellType"
snscc_filtered_regionstats <- RegionStats(snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)
Idents(snscc_filtered_regionstats) <- "CellType_sample"
levels(snscc_filtered_regionstats) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")

#LinkstoPeaks <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA")
#saveRDS(object = LinkstoPeaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_regionstats_CellType_sample_LinkstoPeaks.rds")
#LinkstoPeaks_links <- Links(object = LinkstoPeaks)
#saveRDS(object = LinkstoPeaks_links, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_regionstats_CellType_sample_LinkstoPeaks_peak_rna_trasfer_full_avg_link.rds")
#write.table(LinkstoPeaks_links, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_LinkstoPeaks_links.txt", sep = "\t")

# version 2
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(Normal)","Fibroblasts(Normal)","Endothelial cells(Normal)","Myeloid cells(Normal)","T cells(Normal)","Plasma cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")
snscc_filtered_regionstats <- RegionStats(snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

LinkstoPeaks <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA")
saveRDS(object = LinkstoPeaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_regionstats_CellType_sample_LinkstoPeaks_2.rds")
LinkstoPeaks_links <- Links(object = LinkstoPeaks)
saveRDS(object = LinkstoPeaks_links, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_regionstats_CellType_sample_LinkstoPeaks_peak_rna_trasfer_full_avg_link_2.rds")
write.table(LinkstoPeaks_links, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_LinkstoPeaks_links_2.txt", sep = "\t")

LinkstoPeaks_links_score_cutoff_0 <- LinkstoPeaks_links[LinkstoPeaks_links$score > 0 & LinkstoPeaks_links$pvalue < 0.05, ]
write.table(LinkstoPeaks_links_score_cutoff_0, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15LinkstoPeaks_links_score_cutoff_0_2.txt", sep = "\t")

# heatmap
library(pheatmap)
library(RColorBrewer)
library(dplyr)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")
LinkstoPeaks_links <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_LinkstoPeaks_links.txt")
LinkstoPeaks_links.average <- AverageExpression(snscc_filtered, features = LinkstoPeaks_links$peaks, assays = "peaks")
pheatmap::pheatmap(LinkstoPeaks_links.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, fontsize_row = 5, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058_heatmap_7*8.pdf", plot = figure, height = 8, width = 7, units = "in")

LinkstoPeaks_links.average <- AverageExpression(snscc_filtered, features = LinkstoPeaks_links$gene, assays = "RNA")
pheatmap::pheatmap(LinkstoPeaks_links.average[["RNA"]], scale = 'row', cluster_rows = T, cluster_cols = F, show_rownames = F, fontsize_row = 5, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DEPs_logFC_058_minpct_001_ovelap_DEGs_logFC_058_heatmap_7*8.pdf", plot = figure, height = 8, width = 7, units = "in")

gene_name <- "ENO1"
snscc_filtered_motif <- LinkPeaks(snscc_filtered_motif, peak.assay = "peaks", expression.assay = "RNA", min.cells = 3, score_cutoff = 0.05, genes.use = gene_name, pvalue_cutoff = 0.05)
linked_gene <- GetLinkedPeaks(snscc_filtered_motif, gene_name, min.abs.score = 0, assay = "peaks") # min.abs.score = 0.1

library(ggplot2)
# link peaks to genes (Epi normal)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("PIGR")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F, annotation = T, links = F, extend.upstream = -10000, extend.downstream = 10000, idents = c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_PIGR_Epi_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")
# link peaks to genes (Epi tumor)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("KRT17")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F, annotation = T, links = F, extend.upstream = 1000, extend.downstream = 5000, idents = c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_KRT17_Epi_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")

# link peaks to genes (Fibro normal)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("MYH11")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F, annotation = T, links = F, extend.upstream = -150000, extend.downstream = 5000, idents = c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_MYH11_Fibro_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")
# link peaks to genes (Fibro tumor)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("FAP")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F , annotation = T, links = F, extend.upstream = 2000, extend.downstream = 30000, idents = c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_FAP_Fibro_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")

# link peaks to genes (Endo normal)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("CCL14")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F, annotation = T, links = F, extend.upstream = 1000, extend.downstream = 3000, idents = c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_CCL14_Endo_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")
# link peaks to genes (Endo tumor)
color.track <- c("#D05353","#D05353","#EA9014","#EA9014","#FCCF55","#FCCF55")
test_gene <- c("ESM1")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = F, annotation = T, links = F, extend.upstream = 15000, extend.downstream = 1000, idents = c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)")) & scale_fill_manual(values = color.track) & theme(axis.title.y = element_text(size = 9)) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_ESM1_Endo_7*4.pdf", plot = figure, width = 7, height = 4, units = "in")





# link peaks to genes
test_gene <- c("ADM")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 1000, extend.downstream = 5000, group.by = "sample") & scale_fill_manual(values = color.sample)
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 1000, extend.downstream = 5000, group.by = "CellType_sample", ymax = 40) #& scale_fill_manual(color.snscc.filtered.CellType.sample)


test_gene <- c("KRT6A")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 1000, extend.downstream = 5000, group.by = "sample") & scale_fill_manual(values = color.sample)
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 1000, extend.downstream = 5000, group.by = "CellType_sample", ymax = 40) #& scale_fill_manual(color.snscc.filtered.CellType.sample)

# link peaks to genes
test_gene <- c("ADM")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 5000, extend.downstream = 1500, group.by = "sample") & scale_fill_manual(values = color.sample) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_ADM_sample_9*6.pdf", plot = figure, width = 9, height = 6, units = "in")

test_gene <- c("KRT6A")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 1000, extend.downstream = 5000, group.by = "sample") & scale_fill_manual(values = color.sample) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_res_15_CellType_Track_KRT6A_sample_9*6.pdf", plot = figure, width = 9, height = 6, units = "in")

library(scCustomize)
DefaultAssay(snscc_filtered) <- "RNA"
FeaturePlot_scCustom(snscc_filtered, features = test_gene, reduction = "umap.wnn", colors_use = rev(hcl.colors(5,"Reds")), pt.size = 0.01, split.by = "dataset")

ranges.show <- StringToGRanges("chr11-18407401-18407766")
ranges.show$color <- "orange"
CoveragePlot(object = snscc_filtered_regionstats, region = "chr11-18407401-18407766", region.highlight = ranges.show, peaks = T, annotation = T, links = T, extend.upstream = 5000, extend.downstream = 5000, features = "LDHA") & scale_fill_manual(values = color.snscc.filtered.CellType)

# link peaks to genes (CellType)
# c("CDH1","PDGFRB","CDH5","ITGAX","CD2","FCRL5")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CDH1")
CoveragePlot(object = snscc_filtered_regionstats, region = "CDH1", extend.upstream = 6000, extend.downstream = -90000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_CDH1_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "PDGFRB")
CoveragePlot(object = snscc_filtered_regionstats, region = "PDGFRB", extend.upstream = -36000, extend.downstream = 8000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_PDGFRB_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "DCN")
CoveragePlot(object = snscc_filtered_regionstats, region = "DCN", extend.upstream = -35000, extend.downstream = 8000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_DCN_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CD34")
CoveragePlot(object = snscc_filtered_regionstats, region = "CD34", extend.upstream = -5000, extend.downstream = 10000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_CD34_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "VWF")
CoveragePlot(object = snscc_filtered_regionstats, region = "VWF", extend.upstream = -150000, extend.downstream = 20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_VWF_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CDH5")
CoveragePlot(object = snscc_filtered_regionstats, region = "CDH5", extend.upstream = 10000, extend.downstream = -30000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_CDH5_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "ITGAX")
CoveragePlot(object = snscc_filtered_regionstats, region = "ITGAX", extend.upstream = 6000, extend.downstream = -20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_ITGAX_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CD2")
CoveragePlot(object = snscc_filtered_regionstats, region = "CD2", extend.upstream = 18000, extend.downstream = -10000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_CD2_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "CD3D")
CoveragePlot(object = snscc_filtered_regionstats, region = "CD3D", extend.upstream = 800, extend.downstream = 2500, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_CD3D_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")

snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = "FCRL5")
CoveragePlot(object = snscc_filtered_regionstats, region = "FCRL5", extend.upstream = 10000, extend.downstream = -20000, peaks = T, annotation = T, links = F) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_filtered_res_08_CellType_Track_FCRL5_4*6.pdf", plot = figure, height = 6, width = 4, units = "in")


test_gene <- c("LAMB3")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, peaks = T, annotation = T, links = T, extend.upstream = -25000, extend.downstream = 8000, features = test_gene) & scale_fill_manual(values = color.snscc.filtered.CellType) -> figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/Track_LAMB3_6*4.5.pdf", plot = figure, height = 4.5, width = 6, units = "in")
FeaturePlot_scCustom(snscc_filtered, features = c(test_gene), reduction = "umap.wnn", colors_use = c("#DDDDDD","#DDDDDD","#FFAA93","#FF755F","#EB2C31","#CC1A36","#AC0535","#6D0026"), pt.size = 0.01) -> figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/LAMB3&DSG3_featureplot_8.5*4.pdf", plot = figure, height = 4, width = 8.5, units = "in")


###
## Motif analysis
library(Signac)
library(Seurat)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

DefaultAssay(snscc_filtered) <- "peaks"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))

# add motif information
snscc_filtered_motif <- AddMotifs(object = snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm, assay = "peaks", verbose = TRUE)

# Scan the DNA sequence of each peak for the presence of each motif
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.matrix <- CreateMotifMatrix(features = StringToGRanges(rownames(snscc_filtered_motif), sep = c(":", "-")), pwm = pfm, genome = 'BSgenome.Hsapiens.UCSC.hg38', sep = c(":", "-"))
motif <- CreateMotifObject(data = motif.matrix, pwm = pfm)

# Computing motif activities
snscc_filtered_motif <- RunChromVAR(object = snscc_filtered_motif, genome = BSgenome.Hsapiens.UCSC.hg38)

saveRDS(object = snscc_filtered_motif, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_res_15_motif_CellType_sample.rds")

snscc_filtered_motif <- readRDS("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/RDS_files/snscc_filtered_res_15_motif_CellType_sample.rds")

## Motif analysis - all open peaks
DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered_open_peaks <- AccessiblePeaks(snscc_filtered)
#write.table(snscc_filtered_open_peaks, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_all_open_peaks.txt", sep = "\t")
snscc_filtered_open_peaks <- read.table("/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_all_open_peaks.txt", sep = "\t")
da_peaks <- snscc_filtered_open_peaks$x
peaks_annotation <- ClosestFeature(snscc_filtered, da_peaks)
head(peaks_annotation)
#write.table(peaks_annotation, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_all_open_peaks_annotation.txt", sep = '\t')

# DA motifs (CellType)
DefaultAssay(snscc_filtered_motif) <- 'chromvar'
Idents(snscc_filtered_motif) <- "CellType"
library(pheatmap)
library(RColorBrewer)
library(dplyr)
differential.activity <- FindAllMarkers(object = snscc_filtered_motif, only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 1)
differential.activity.cut <- differential.activity[differential.activity$p_val_adj < 0.05 & differential.activity$avg_log2FC > 1, ]
differential.activity.cut <- differential.activity[differential.activity$p_val_adj < 0.05, ]
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = differential.activity.cut$gene, assays = "chromvar")

motifLookup <- differential.activity.cut$gene
motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
differential.activity.cut$gene_symbol <- motifNames
write.table(differential.activity.cut, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_motif_CellType_DaTFs_log2FC_058_rowmeans.txt", sep = '\t')

pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, annotation_col = anno_col, annotation_colors = anno_colors, labels_row = differential.activity.cut$gene_symbol, fontsize_row = 5, color = colorRampPalette(rev(brewer.pal(n=7,name="RdYlGn")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_DaTFs_logFC_058_7*40.pdf", plot = figure, height = 40, width = 7, units = "in")

# DA motifs (CellType_sample)
Idents(snscc_filtered_motif) <- "CellType_sample"
#levels(snscc_filtered_motif) <- c("Plasma cells(Tumor)","Plasma cells(Normal)","T cells(Tumor)","T cells(Normal)","Myeloid cells(Tumor)","Myeloid cells(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Epithelial cells(Tumor)","Epithelial cells(Normal)")
levels(snscc_filtered_motif) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")
DefaultAssay(snscc_filtered_motif) <- "chromvar"
library(pheatmap)
library(RColorBrewer)
library(dplyr)
differential.activity <- FindAllMarkers(object = snscc_filtered_motif, only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.58)
differential.activity.cut <- differential.activity[differential.activity$p_val_adj < 0.05, ]
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = differential.activity.cut$gene, assays = "chromvar")

motifLookup <- differential.activity.cut$gene
motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
differential.activity.cut$gene_symbol <- motifNames
write.table(differential.activity.cut, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_motif_CellType_sample_DaTFs_logFC_1_rowMeans_reverse.txt", sep = '\t')
write.table(differential.activity.cut, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_motif_CellType_sample_DaTFs_logFC_058_rowMeans_reverse.txt", sep = '\t')

anno_col <- snscc_filtered@meta.data %>% data.frame() %>% select(sample,CellType_sample) %>% unique()
anno_col
rownames(anno_col) <- anno_col$CellType_sample
anno_col
anno_colors <- list(CellType_sample = c("Epithelial cells(Tumor)" = "#D05353","Fibroblasts(Tumor)" = "#EA9014","Endothelial cells(Tumor)" = "#FCCF55","Myeloid cells(Tumor)" = "#76AB62","T cells(Tumor)" = "#4d648d","Plasma cells(Tumor)"="#9F8CB3","Epithelial cells(Normal)" = "#D05353","Fibroblasts(Normal)" = "#EA9014","Endothelial cells(Normal)" = "#FCCF55","Myeloid cells(Normal)" = "#76AB62","T cells(Normal)" = "#4d648d","Plasma cells(Normal)"="#9F8CB3"), sample = c("Normal"="#4A6274", "Tumor"="#E2725A"))
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_sample_DaTF_log2FC_1_reverse_4*8.pdf", plot = figure, height = 8, width = 4, units = "in")
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = T, cutree_rows = 12, cluster_cols = F, show_rownames = F, annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100)) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_sample_DaTF_log2FC_1_reverse_cut_tree_12_6*8.pdf", plot = figure, height = 8, width = 6, units = "in")
heatmap_cluster <- cutree(figure$tree_row,12)
heatmap_cluster_dataframe <- data.frame(heatmap_cluster)
#write.table(heatmap_cluster_dataframe, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_CellType_sample_DaTF_log2FC_1_cutree_cluster.txt", sep = "\t")

# DA motifs (sample)
Idents(snscc_filtered_motif) <- "sample"
DefaultAssay(snscc_filtered_motif) <- "chromvar"
library(pheatmap)
library(RColorBrewer)
library(dplyr)
differential.activity <- FindMarkers(object = snscc_filtered_motif, only.pos = F, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 1, ident.1 = "Tumor", ident.2 = "Normal")
differential.activity.cut <- differential.activity[differential.activity$p_val_adj < 0.05, ]
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = rownames(differential.activity.cut), assays = "chromvar")

motifLookup <- rownames(differential.activity.cut)
motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
differential.activity.cut$gene_symbol <- motifNames
write.table(differential.activity.cut, file = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Documents/snscc_filtered_res_15_motif_sample_DaTFs_FindMarkers_log2FC_1_rowMeans.txt", sep = '\t')

Idents(snscc_filtered_motif) <- "CellType_sample"
levels(snscc_filtered_motif) <- c("Epithelial cells(Normal)","Epithelial cells(Tumor)","Fibroblasts(Normal)","Fibroblasts(Tumor)","Endothelial cells(Normal)","Endothelial cells(Tumor)","Myeloid cells(Normal)","Myeloid cells(Tumor)","T cells(Normal)","T cells(Tumor)","Plasma cells(Normal)","Plasma cells(Tumor)")
motif_list <- c("NFIA","NFIB","NFIX","SOX4","SOX10","FOSL1","FOSL2","JUNB","HIF1A","NR3C1","NR3C2","SP2","SP8","NFATC3","NFATC4","RUNX2","TWIST1","GATA2","GATA3","PPARA::RXRA","RXRB","ETV1","ETV4","SPI1","SPIB","IRF1","USF1","USF2","ETS1","ETS2","IRF4")
motif.name <- ConvertMotifID(snscc_filtered_motif, name = motif_list, assay = "peaks")
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = motif.name, assays = "chromvar")
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", labels_row = motif_list, fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=11,name="PiYG")))(100)) -> figure #PiYG #BrBG
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_sample_DaTF_log2FC_058_sorting_heatmap_4*8.pdf", plot = figure, height = 8, width = 4, units = "in")
#library(ArchR)
#pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", labels_row = motif_list, fontsize_row = 8, color = paletteContinuous("solarExtra"))
DefaultAssay(snscc_filtered) <- "RNA"
library(scCustomize)
library(RColorBrewer)
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = motif_list, assays = "RNA")
pheatmap::pheatmap(differential.activity.cut.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100))
FeaturePlot_scCustom(snscc_filtered, features = "GATA4", pt.size = 0.01, split.by = "sample", reduction = "umap.wnn", num_columns = 1, colors_use = (brewer.pal(n=9, name = "Reds")))


motif_list <- c("NFIA","NFIB","NFIX","SOX10","SOX13","TP63","ATF4","CEBPD","FOSL1","FOSL2","JUNB","HIF1A","PRRX1","PRRX2","NFATC2","NFATC4","RUNX2","TWIST1","GATA2","GATA3","PPARA::RXRA","RXRB","NR2F2","ERG","RXRB","SPI1","MITF","USF1","USF2","IRF1","ETS1","IRF4","IRF6")
motif.name <- ConvertMotifID(snscc_filtered_motif, name = motif_list, assay = "peaks")
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = motif.name, assays = "chromvar")
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", labels_row = motif_list, fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=11,name="PiYG")))(100)) -> figure #PiYG #BrBG
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/snscc_filtered_CellType_sample_DaTF_log2FC_058_sorting_heatmap_4*8.pdf", plot = figure, height = 8, width = 4, units = "in")
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = motif_list, assays = "RNA")
pheatmap::pheatmap(differential.activity.cut.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100))


# TP53,TP63,GRHL1,GRHL2,FOSL1,FOSL2/TWIST1,NFATC2,NFATC4,RUNX2,"GLI2","GLI3"/GABPA,ETV1,ETV4,ETS1,ETS2,ERG/SPI1,SPIB,SPIC,"MITF","IRF8"


levels(snscc_filtered_motif) <- c("Epithelial cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Tumor)","Fibroblasts(Normal)","Endothelial cells(Tumor)","Endothelial cells(Normal)","Myeloid cells(Tumor)","Myeloid cells(Normal)","T cells(Tumor)","T cells(Normal)","Plasma cells(Tumor)","Plasma cells(Normal)")

snscc_filtered_motif_tumor <- subset(snscc_filtered_motif, sample == "Tumor")
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif_tumor, features = motif_sorting, assays = "chromvar")
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", annotation_col = anno_col, annotation_colors = anno_colors, labels_row = differential.activity.cut$motif.name, fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=9,name="RdYlGn")))(100))

snscc_filtered_motif_normal <- subset(snscc_filtered_motif, sample == "Normal")
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif_normal, features = motif_sorting, assays = "chromvar")
pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, border_color = "black", annotation_col = anno_col, annotation_colors = anno_colors, labels_row = differential.activity.cut$motif.name, fontsize_row = 8, color = colorRampPalette(rev(brewer.pal(n=9,name="RdYlGn")))(100))

#ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/snscc_filtered_CellType_enriched_motif_custom_reverse_4*8.pdf", plot = figure, height = 8, width = 4, units = "in")
#levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells")

MotifPlot(snscc_filtered_motif, motifs = c("MA0525.2","MA0647.1","MA0478.1"), assay = "peaks") -> a1
MotifPlot(snscc_filtered_motif, motifs = c("MA1123.2","MA0474.2","MA0080.5"), assay = "peaks") -> a2
a3 <- ggarrange(a1,a2,nrow=1)
a3

# look at the activity and gene expression of motif
test_gene <- c("FOSL1")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", pt.size = 0.01, colors_use = c("#002153","#00366C","#386D9E","#8CABC7","#FFEE99","#D07C42","#C7522B","#b21704","#611300"), na_cutoff = NULL) -> p1
p1
#c("#002153","#00366C","#0072B4","#79ABE2","#FFFFBF","#D07C42","#C7522B","#b21704","#611300") #A9BA9D
#paletteContinuous("solarExtra",n=8)
#Plot_Density_Custom(snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", pt.size = 0.01, custom_palette = c("#002153","#00366C","#386D9E","#8CABC7","#FFEE99","#D07C42","#C7522B","#b21704","#611300"))
DefaultAssay(snscc_filtered_motif) <- "RNA"
FeaturePlot_scCustom(snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9, "OrRd"))) -> p2
p2
#Plot_Density_Custom(snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.01, custom_palette = c("#002153","#00366C","#386D9E","#8CABC7","#FFEE99","#D07C42","#C7522B","#b21704","#611300"))
library(ggpubr)
ggarrange(p1,p2,nrow=2)
ggarrange(p1,p2,nrow=2) -> figure
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/FOSL1_featureplot_motif_with_gene_exp_4.5*8.pdf", plot = figure, width = 4.5, height = 8, units = "in")
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/NFIX_featureplot_motif_with_gene_exp_4.5*8.pdf", plot = figure, width = 4.5, height = 8, units = "in")
#ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/RUNX2_featureplot_motif_with_gene_exp_4.5*8.pdf", plot = figure, width = 4.5, height = 8, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/Figures/RUNX2_featureplot_motif_with_gene_exp_4.5*8.pdf", plot = figure, width = 4.5, height = 8, units = "in")

# look at the activity and gene expression of motif
Idents(snscc_filtered_motif) <- "CellType_sample"
levels(snscc_filtered_motif) <- c("Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)","Epithelial cells(Normal)","Fibroblasts(Normal)","Endothelial cells(Normal)","Myeloid cells(Normal)","T cells(Normal)","Plasma cells(Normal)")
test_gene <- c("HIF1A")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
color.snscc.filtered.CellType.sample <- c("#D05353", "#EA9014", "#FCCF55", "#76AB62", "#4d648d", "#9F8CB3", "#D05353", "#EA9014", "#FCCF55", "#76AB62", "#4d648d", "#9F8CB3")
VlnPlot(snscc_filtered_motif, features = motif.name, pt.size = 0, cols = color.snscc.filtered.CellType.sample) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() -> p1
DefaultAssay(snscc_filtered_motif) <- "RNA"
VlnPlot(snscc_filtered_motif, features = test_gene, pt.size = 0, cols = color.snscc.filtered.CellType.sample) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() -> p2
library(ggpubr)
ggarrange(p1,p2,nrow=1)
p1
p2
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/HIF1A_featureplot_motif_10*4.pdf", plot = p1, width = 10, height = 4, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/HIF1A_featureplot_gene_exp_9*4.pdf", plot = p2, width = 9, height = 4, units = "in")

# look at the activity and gene expression of motif
Idents(snscc_filtered_motif) <- "sample"
test_gene <- c("TP63")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
color.sample <- c("#4A6274", "#E2725A")
VlnPlot(snscc_filtered_motif, features = motif.name, pt.size = 0, cols = color.sample) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() -> p1
DefaultAssay(snscc_filtered_motif) <- "RNA"
VlnPlot(snscc_filtered_motif, features = test_gene, pt.size = 0, cols = color.sample) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.1) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() -> p2
library(ggpubr)
ggarrange(p1,p2,nrow=1) -> figure
figure
ggsave(filename = "/home/iphd2/Desktop/SNSCC_Normal_multiome_2023/SoupX/Figures/HIF1A_featureplot_motif_10*4.pdf", plot = figure, width = 10, height = 4, units = "in")

# look at the activity and gene expression of motif
# TP63, GRHL2, TWIST1, NFATC2, ERG, SPI1
test_gene <- c("FOSL1","TP63","RUNX2","TWIST1","ERG","SPI1")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
#motif.name <- c("MA1105.2", "MA0525.2", "MA0511.2", "MA0633.1", "MA0474.2","MA1484.1","MA0080.5","MA1420.1")
FeaturePlot_scCustom(snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", pt.size = 0.01, num_columns = 6, colors_use = c("#002153","#00366C","#0072B4","#79ABE2","#FFFFBF","#D07C42","#C7522B","#b21704","#611300"), na_cutoff = NULL) -> p1
DefaultAssay(snscc_filtered_motif) <- "RNA"
FeaturePlot_scCustom(snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.01, num_columns = 6, colors_use = (brewer.pal(n=9, "Reds"))) -> p2
library(ggpubr)
ggarrange(p1,p2,nrow=2)
ggarrange(p1,p2,nrow=2) -> figure
#a3 -> figure
ggsave(filename = "/home/iphd2/Desktop/snscc_multiome_CA5_CA8/Figure_save_update/Enriched_motif_with_gene_exp_custom_26*8.pdf", plot = figure, height = 8, width = 26, units = "in")


# look at the activity and gene expression of JUNB
library(scCustomize)
test_gene <- "MEIS3"
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", colors_use = brewer.pal(n=9,"Reds"), pt.size = 0.01)
#VlnPlot(snscc_filtered_motif, features = motif.name, pt.size = 0.1, split.by = "dataset")
DefaultAssay(snscc_filtered_motif) <- "RNA"
FeaturePlot_scCustom(snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.01, colors_use = (brewer.pal(n=9,"Reds")))
#library(ggpubr)
#ggarrange(p1,p2,nrow=1)
VlnPlot(snscc_filtered_motif, features = "FOXA2", pt.size = 0.1, split.by = "dataset")

# look at the activity and gene expression of JUNB
test_gene <- "IKZF1"
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", pt.size = 0.01) -> p1
DefaultAssay(snscc_filtered_motif) <- "RNA"
FeaturePlot_scCustom(snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.01) -> p2
library(ggpubr)
ggarrange(p1,p2,nrow=1)

# look at the activity of JUNB
DefaultAssay(snscc_filtered_motif) <- "peaks"
p1 <- DimPlot(snscc_filtered_motif, label = T, reduction = "umap.wnn", pt.size = 0.1) + NoLegend()
motif.name <- ConvertMotifID(snscc_filtered_motif, name = 'TP63')
DefaultAssay(snscc_filtered_motif) <- "chromvar"
p2 <- FeaturePlot(object = snscc_filtered_motif, features = motif.name, min.cutoff = 'q10', max.cutoff = 'q90', reduction = "umap.wnn", pt.size = 0.1, split.by = "dataset") #can replace cutoff option to min.cutoff = 0
p1+p2

## Footprint
# gather the footprinting information for sets of motifs
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif_name <- c("JUNB")
snscc_filtered_motif <- Footprint(
  object = snscc_filtered_motif,
  motif.name = motif_name,
  genome = BSgenome.Hsapiens.UCSC.hg38)

# plot the footprint data for each group of cells
PlotFootprint(snscc_filtered_motif, features = motif_name) & scale_fill_manual(values = color.snscc.filtered.CellType)

#############################################################################3

# look at the activity and gene expression of JUNB
#test_gene <- "JUNB"
#DefaultAssay(snscc_filtered_motif) <- "peaks"
#motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
#DefaultAssay(snscc_filtered_motif) <- "chromvar"
#Plot_Density_Custom(seurat_object = snscc_filtered_motif, features = motif.name, reduction = "umap.wnn", pt.size = 0.1, custom_palette = viridis_light_high) -> p1
#DefaultAssay(snscc_filtered_motif) <- "RNA"
#Plot_Density_Custom(seurat_object = snscc_filtered_motif, features = test_gene, reduction = "umap.wnn", pt.size = 0.1, custom_palette = viridis_light_high) -> p2
#library(ggpubr)
#ggarrange(p1,p2,nrow=1)

# look at the activity of JUNB
#DefaultAssay(snscc_filtered_motif) <- "peaks"
#p1 <- DimPlot(snscc_filtered_motif, label = T, reduction = "umap.wnn", pt.size = 0.1) + NoLegend()
#motif.name <- ConvertMotifID(snscc_filtered_motif, name = 'TP63')
#DefaultAssay(snscc_filtered_motif) <- "chromvar"
#p2 <- FeaturePlot(object = snscc_filtered_motif, features = motif.name, min.cutoff = 'q10', max.cutoff = 'q90', reduction = "umap.wnn", pt.size = 0.1, split.by = "dataset") #can replace cutoff option to min.cutoff = 0
#p1+p2

