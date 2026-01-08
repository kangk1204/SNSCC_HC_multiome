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
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("KRT6A"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS2D_KRT6A
FigureS2D_KRT6A
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_KRT6A_4*9.pdf", plot = FigureS2D_KRT6A, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NDRG1"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS2D_NDRG1
FigureS2D_NDRG1
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_NDRG1_4*9.pdf", plot = FigureS2D_NDRG1, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FAP"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS2D_FAP
FigureS2D_FAP
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_FAP_4*9.pdf", plot = FigureS2D_FAP, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("ANGPT2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS2D_ANGPT2
FigureS2D_ANGPT2
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_ANGPT2_4*9.pdf", plot = FigureS2D_ANGPT2, height = 4, width = 9, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("IL2RB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS2D_IL2RB
FigureS2D_IL2RB
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_IL2RB_4*9.pdf", plot = FigureS2D_IL2RB, height = 4, width = 9, units = "in")



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