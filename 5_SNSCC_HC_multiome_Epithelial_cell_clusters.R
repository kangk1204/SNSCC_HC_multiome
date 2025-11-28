library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Upload RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")



# Subset Epithelial cell cluster
Epi <- subset(snscc_filtered, CellType == "Epithelial cells")
DimPlot(Epi, label = TRUE, reduction = "umap.wnn", group.by = "dataset", pt.size = 0.01)



# Gene expression data processing
library(harmony)
DefaultAssay(Epi) <- "RNA"
Epi <- FindVariableFeatures(Epi, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
Epi <- ScaleData(Epi, features=rownames(Epi))
Epi <- RunPCA(Epi)
DepthCor(Epi, n = 50, reduction = "harmony_pca")
ElbowPlot(Epi, ndims = 50, reduction = "harmony_pca")
Epi <- RunUMAP(Epi, dims = 1:30, reduction = "harmony_pca", reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DimPlot(object = Epi, label = TRUE, reduction = "umap.rna", group.by = "dataset")



# DNA accessibility data processing
DefaultAssay(Epi) <- "peaks"
Epi <- FindTopFeatures(Epi, min.cutoff = 5)
Epi <- RunTFIDF(Epi)
Epi <- RunSVD(Epi)
DepthCor(Epi, n = 50, reduction = "harmony_lsi")
ElbowPlot(Epi, ndims = 50, reduction = "harmony_lsi")
Epi <- RunUMAP(Epi, reduction = "harmony_lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(object = Epi, label = TRUE, reduction = "umap.atac", group.by = "dataset")



# Joint UMAP visualization
# build a joint neighbor graph using both assays
Epi <- FindMultiModalNeighbors(
  object = Epi,
  reduction.list = list("harmony_pca", "harmony_lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = c("RNA.weight","peaks.weight"),
  verbose = TRUE) 



# build a joint UMAP visualization
Epi <- RunUMAP(
  object = Epi,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "umap.wnn",
  reduction.key = "wnnUMAP_",
  verbose = TRUE)



# Clustering
Epi <- FindClusters(Epi, resolution = 0.45, algorithm = 3, graph.name = "wsnn")
DimPlot(Epi, reduction = "umap.wnn", pt.size = 0.05, label = T)



# Cell type labeling
Epi_rename <- RenameIdents(Epi,"0"="Serous","1"="TC1","2"="Serous","3"="Serous","4"="TC2","5"="TC3","6"="TC4","7"="Basal","8"="Goblet","9"="Mucous","10"="Ciliated","11"="Serous","12"="TC5","13"="MEC","14"="Tuft")
Epi_rename$CellType_epi <- Idents(Epi_rename)
Idents(Epi_rename) <- "CellType_epi"
levels(Epi_rename) <- c("Basal","Ciliated","Goblet","MEC","Mucous","Serous","Tuft","TC1","TC2","TC3","TC4","TC5")
DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, label = T)
DimPlot(Epi_rename, reduction = "umap.rna", pt.size = 0.05, label = T)
DimPlot(Epi_rename, reduction = "umap.atac", pt.size = 0.05, label = T)

Epi_rename$CellType_epi_sample <- paste(Epi_rename$CellType_epi,"_",Epi_rename$sample)
DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, label = T, group.by = "CellType_epi_sample")

Epi_rename$CellType_epi_dataset <- paste(Epi_rename$CellType_epi,"_",Epi_rename$dataset)
DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, label = T, group.by = "CellType_epi_dataset")

saveRDS(Epi_rename, file = "./RDS_files/Epi_rename_res_045.rds")



# UMAP plots
Epi_rename <- readRDS("./RDS_files/Epi_rename_res_045.rds")

color.Epi.rename <- c("#86B6CB","#277AA5","#194F6F","#C6B1CF","#89608E","#757C98","#B3B5AA","#78010A","#E10E07","#EC796E","#E3989C","#FBC8BB")
color.sample <- c("#F7C475","#7D5008")

DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, group.by = "seurat_clusters") -> FigureS5A
FigureS5A
ggsave(filename = "./Figures/Epi_rename_res_045_umap_wnn_4.5*4.pdf", plot = FigureS5A, width = 4.5, height = 4, units = "in")

DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, group.by = "dataset") -> FigureS5B
FigureS5B
ggsave(filename = "./Figures/Epi_rename_res_045_CellType_epi_umap_wnn_groupby_dataset_4.8*4.pdf", plot = FigureS5B, width = 4.8, height = 4, units = "in")

DimPlot(Epi_rename, reduction = "umap.wnn", pt.size = 0.05, split.by = "sample", cols = color.Epi.rename) -> Figure5A
Figure5A
ggsave(filename = "./Figures/Epi_rename_res_045_CellType_epi_umap_wnn_splitby_sample_8*4.pdf", plot = Figure5A, width = 8, height = 4, units = "in")



# Cell proportion
library(dittoSeq)
dittoBarPlot(Epi_rename, var = "CellType_epi", group.by = "dataset", x.labels.rotate = F, color.panel = color.Epi.rename, var.labels.reorder = c(1,2,3,4,5,6,12,7,8,9,10,11), x.reorder = c(3,4,1,2)) -> Figure5B
Figure5B
ggsave(filename = "./Figures/Epi_rename_res_045_CellType_epi_dataset_proportion_3*4.pdf", plot = Figure5B, width = 3, height = 4, units = "in")

dittoBarPlot(Epi_rename, var = "sample", group.by = "CellType_epi", x.labels.rotate = T, color.panel = color.sample, x.reorder = c(1,2,3,4,5,6,12,7,8,9,10,11)) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_continuous(position = "right") -> FigureS5D
FigureS5D
ggsave(filename = "./Figures/Epi_rename_res_045_CellType_epi_proportion_6*3.pdf", plot = FigureS5D, width = 6, height = 3, units = "in")



# Marker gene expression
DefaultAssay(Epi_rename) <- "RNA"
Idents(Epi_rename) <- "seurat_clusters"
DotPlot(Epi_rename, features = c("KRT15","TP63","FOXJ1","PIFO","MUC5AC","SCGB1A1","ACTA2","MYH11","MUC5B","BPIFB2","LYZ","LTF","POU2F3","TRPM5"), cols = c("#DDDDDD","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> FigureS5C
FigureS5C
ggsave(filename = "./Figures/Epi_rename_res_045_15clusters_normal_gene_exp_dotplot_5.5*4.pdf", plot = FigureS5C, width = 5.5, height = 4, units = "in")

Idents(Epi_rename) <- "CellType_epi"
levels(Epi_rename) <- c("Basal","Ciliated","Goblet","MEC","Mucous","Serous","Tuft","TC1","TC2","TC3","TC4","TC5")
DotPlot(Epi_rename, features = c("AZGP1","PROM1","KRT6A","SFN")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + scale_colour_gradientn(colours = c("yellow","#d896ff","#be29ec","#660066")) -> Figure5C
Figure5C
ggsave(filename = "./Figures/Epi_rename_res_045_tumor_gene_exp_dotplot_3.5*3.5.pdf", plot = Figure5C, width = 3.5, height = 3.5, units = "in")



# AddModuleScore
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scCustomize)

AddModuleScore_gene_list <- read.csv("./Data/AddModuleScore_2024.csv", header = T, sep = ",")

Module <- list(AddModuleScore_gene_list$bulk_down_DEGs_top100)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Down")
Module <- list(AddModuleScore_gene_list$bulk_up_DEGs_top100)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Up")

VlnPlot(Epi_rename, features = "Down1", pt.size = 0, cols = color.Epi.rename) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() & theme(text = element_text(size = 20)) -> FigureS5E_down
FigureS5E_down
ggsave(filename = "./Figures/Epi_rename_res_045_bulk_Down_DEGs_Top_100_module_score_vlnplot_7*4.pdf", plot = FigureS5E_down, height = 4, width = 7, units = "in")

VlnPlot(Epi_rename, features = "Up1", pt.size = 0, cols = color.Epi.rename) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.15) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() + RotatedAxis() & theme(text = element_text(size = 20)) -> FigureS5E_up
FigureS5E_up
ggsave(filename = "./Figures/Epi_rename_res_045_bulk_Up_DEGs_Top_100_module_score_vlnplot_7*4.pdf", plot = FigureS5E_up, height = 4, width = 7, units = "in")

AddModuleScore_gene_list <- read.csv("./Data/AddModuleScore_HNSCC.csv", header = T, sep = ",")

Module <- list(AddModuleScore_gene_list$Cell_cycle)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Cell_cycle")
Module <- list(AddModuleScore_gene_list$pEMT)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "pEMT")
Module <- list(AddModuleScore_gene_list$Epi_diff1)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Epi_diff1")
Module <- list(AddModuleScore_gene_list$Epi_diff2)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Epi_diff2")
Module <- list(AddModuleScore_gene_list$Epi_diff)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Epi_diff")
Module <- list(AddModuleScore_gene_list$Stress)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Stress")
Module <- list(AddModuleScore_gene_list$Hypoxia)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Hypoxia")

DotPlot(Epi_rename, features = c("Hypoxia1","Stress1","Cell_cycle1","Epi_diff11","Epi_diff21","pEMT1"), cols = "RdBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> FigureS5J
FigureS5J
ggsave(filename = "./Figures/HNSCC_gene_list_modulescore_4*4.pdf", plot = FigureS5J, width = 4, height = 4, units = "in")



# DEGs (Tumor cell clusters)
Epi_rename <- RenameIdents(Epi_rename, "TC1"="TC","TC2"="TC","TC3"="TC","TC4"="TC","TC5"="TC")
Epi_rename$TC <- Idents(Epi_rename)
Idents(Epi_rename) <- "CellType_epi"
TC <- subset(Epi_rename, TC == "TC")
Idents(TC) <- "CellType_epi"
TC$CellType_epi <- Idents(TC)

color.TC <- c("#78010A","#E10E07","#EC796E","#E3989C","#FBC8BB")

# DEGs heatmap
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
DefaultAssay(TC) <- "RNA"
levels(TC) <- c("TC1","TC2","TC3","TC4","TC5")
TC.all.genes <- FindAllMarkers(TC, min.pct = 0.25, only.pos = T, logfc.threshold = 0.58)
TC.all.da.genes <- TC.all.genes[TC.all.genes$p_val_adj < 0.05 & TC.all.genes$avg_log2FC > 0.58, ]
TC.da.genes.average <- AverageExpression(TC, features = TC.all.da.genes$gene, assays = "RNA")
pheatmap::pheatmap(TC.da.genes.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)) -> figure
ggsave(filename = "./Figures/Epi_rename_res_045_TC_heatmap_DEGs_log2FC_058_3*6.pdf", plot = figure, width = 3, height = 6, units = "in")
write.table(TC.all.da.genes, file = "./Documents/TC_DEGs_log2FC_058.txt", sep = '\t')

Epi.all.da.genes <- c("ADM","HK2","LDHA","NDRG1","SLC2A1","ATF4","JUN","JUNB","GCLC","PRDX1","AKR1C3","KRT5","KRT17","SFN","CCNB1","MKI67","TOP2A","ACACA","CHD9","NCOA2","NCOA3","CST3","CTSB","CTSD","CXCL17","MUC1","MUC4","FN1","ITGA2","ITGA3","TGFBI","ZEB1")
Epi.da.genes.average <- AverageExpression(TC, features = Epi.all.da.genes, assays = "RNA")
pheatmap::pheatmap(Epi.da.genes.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, fontsize = 7, border_color = "black", color = colorRampPalette(rev(brewer.pal(n=7,name="RdBu")))(100)) -> Figure5G
Figure5G
ggsave(filename = "./Figures/Epi_rename_res_045_TC_custom_heatmap_2*6.pdf", plot = Figure5G, width = 2, height = 6, units = "in")



# Upset Plot
DefaultAssay(TC) <- "RNA"
up_DEGs_AMC <- read.table(file = "./Data/up_DEGs_AMC.txt", sep = '\t', header = T)
TC.all.da.genes <- read.table(file = "./Documents/TC_DEGs_log2FC_058.txt", sep = '\t')
Epi.all.da.genes <- TC.all.da.genes[TC.all.da.genes$gene %in% up_DEGs_AMC$Gene, ]
Epi_DEG_list <- split(Epi.all.da.genes$gene, Epi.all.da.genes$cluster)
list_to_matrix(Epi_DEG_list)
library(ComplexHeatmap)
Epi_DEG_mt <- ComplexHeatmap::make_comb_mat(Epi_DEG_list)
Epi_DEG_mt
ComplexHeatmap::UpSet(Epi_DEG_mt)
ComplexHeatmap:: UpSet(Epi_DEG_mt, set_order = c("TC1","TC2","TC3","TC4","TC5"), row_names_gp=gpar(fontsize=10), right_annotation = upset_right_annotation(Epi_DEG_mt, gp=gpar(fill=color.TC), add_numbers=T, width = unit(2,"cm"), annotation_name_gp=gpar(fontsize=10)), top_annotation = upset_top_annotation(Epi_DEG_mt, add_numbers=T, annotation_name_gp=gpar(fontsize=10))) -> FigureS5H
FigureS5H
ggsave(filename = "./Figures/Epi_rename_res_045_TC_upsetplot_4*3.pdf", plot = FigureS5H, width = 2, height = 6, units = "in")


# Gene expression & track
# link peaks to genes
library(ggplot2)
DefaultAssay(Epi_rename) <- "peaks"
Epi.regionstats <- RegionStats(Epi_rename, genome = BSgenome.Hsapiens.UCSC.hg38)

# TC1
library(ggplot2)
test_gene <- c("NDRG1")
Epi.regionstats <- LinkPeaks(object = Epi.regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = Epi.regionstats, region = test_gene, peaks = T, annotation = T, links = T, extend.upstream = -25000, extend.downstream = 35000, expression.assay = "RNA", features = test_gene, ymax = 25) & scale_fill_manual(values = color.Epi.rename)-> Figure5H_NDRG1
Figure5H_NDRG1
ggsave(filename = "./Figures/Epi_Track_NDRG1_4*6.pdf", plot = Figure5H_NDRG1, width = 4, height = 6, units = "in")

# TC2
test_gene <- c("KRT17")
Epi.regionstats <- LinkPeaks(object = Epi.regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = Epi.regionstats, region = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 8000, extend.downstream = 6000, expression.assay = "RNA", features = test_gene, ymax = 25) & scale_fill_manual(values = color.Epi.rename)-> Figure5H_KRT17
Figure5H_KRT17
ggsave(filename = "./Figures/Epi_Track_KRT17_4*6.pdf", plot = Figure5H_KRT17, width = 4, height = 6, units = "in")

# Gene expression
DefaultAssay(Epi_rename) <- "RNA"

FeaturePlot_scCustom(seurat_object = Epi_rename, features = c("NDRG1"), pt.size = 0.05, reduction = "umap.wnn", colors_use = (brewer.pal(n=9,"OrRd"))) -> Figure5H_NDRG1_gene
Figure5H_NDRG1_gene
ggsave(filename = "./Figures/Epi_rename_res_045_NDRG1_exp_featureplot_4.5*4.pdf", plot = Figure5H_NDRG1_gene, width = 4.5, height = 4, units = "in")

FeaturePlot_scCustom(seurat_object = Epi_rename, features = c("KRT17"), pt.size = 0.05, reduction = "umap.wnn", colors_use = (brewer.pal(n=9,"OrRd"))) -> Figure5H_KRT17_gene
Figure5H_KRT17_gene
ggsave(filename = "./Figures/Epi_rename_res_045_KRT17_exp_featureplot_4.5*4.pdf", plot = Figure5H_KRT17_gene, width = 4.5, height = 4, units = "in")



# GSVA scores
library(GSVA)

setwd("./AMC/")
fpkm_data <- as.matrix(read.table("fpkm_AMC.txt", header = TRUE, sep = "\t", row.names = 1))

gene_list_TC1 <- read.table("TC1.txt", header = TRUE, sep = "\t")
gene_set_TC1 <- list(gene_list_TC1$TC1)

gene_list_TC2 <- read.table("TC2.txt", header = TRUE, sep = "\t")
gene_set_TC2 <- list(gene_list_TC2$TC2)

gene_list_TC3 <- read.table("TC3.txt", header = TRUE, sep = "\t")
gene_set_TC3 <- list(gene_list_TC3$TC3)

gene_list_TC4 <- read.table("TC4.txt", header = TRUE, sep = "\t")
gene_set_TC4 <- list(gene_list_TC4$TC4)

gene_list_TC5 <- read.table("TC5.txt", header = TRUE, sep = "\t")
gene_set_TC5 <- list(gene_list_TC5$TC5)

# Calculating module scores
method <- "gsva"
verbose <- FALSE
gsva_result <- gsva(fpkm_data, gene_set_TC1, method = method, verbose = verbose)

gsva_df <- data.frame(GSVA_score = gsva_result[,])
gsva_df
write.table(gsva_df, file = "TC1_AMC_results.txt", sep = "\t", quote = FALSE, row.names = TRUE) # Figure5I_TC1

gsva_result <- gsva(fpkm_data, gene_set_TC2, method = method, verbose = verbose)
gsva_df <- data.frame(GSVA_score = gsva_result[,])
gsva_df
write.table(gsva_df, file = "TC2_AMC_results.txt", sep = "\t", quote = FALSE, row.names = TRUE) # Figure5I_TC2

gsva_result <- gsva(fpkm_data, gene_set_TC3, method = method, verbose = verbose)
gsva_df <- data.frame(GSVA_score = gsva_result[,])
gsva_df
write.table(gsva_df, file = "TC3_AMC_results.txt", sep = "\t", quote = FALSE, row.names = TRUE) # FigureS5I_TC3

gsva_result <- gsva(fpkm_data, gene_set_TC4, method = method, verbose = verbose)
gsva_df <- data.frame(GSVA_score = gsva_result[,])
gsva_df
write.table(gsva_df, file = "TC4_AMC_results.txt", sep = "\t", quote = FALSE, row.names = TRUE) # FigureS5I_TC4

gsva_result <- gsva(fpkm_data, gene_set_TC5, method = method, verbose = verbose)
gsva_df <- data.frame(GSVA_score = gsva_result[,])
gsva_df
write.table(gsva_df, file = "TC5_AMC_results.txt", sep = "\t", quote = FALSE, row.names = TRUE) # FigureS5I_TC5



## AddModuleScore
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scCustomize)
DefaultAssay(Epi_rename) <- "RNA"

AddModuleScore_gene_list <- AddModuleScore_gene_list <- read.csv("./Data/AddModuleScore_2024.csv", header = T, sep = ",")

# MsigDB hallmark HYPOXIA gene list
Module <- list(AddModuleScore_gene_list$MsigDB_hypoxia)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Hallmark_HYPOXIA")
FeaturePlot_scCustom(Epi_rename, features = "Hallmark_HYPOXIA1", pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#BBBBBB","#CCCCCC","#DDDDDD","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D"), na_cutoff = NULL) -> FigureS6B_Hypoxia
FigureS6B_Hypoxia
ggsave(filename = "./Figures/Epi_rename_res_045_Module_Hallmark_HYPOXIA_wnnUMAP_4.5*4.pdf", plot = FigureS6B_Hypoxia, width = 4.5, height = 4, units = "in")

VlnPlot(Epi_rename, features = "Hallmark_HYPOXIA1", pt.size = 0, cols = color.Epi.rename) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.2) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() -> FigureS6B_Hypoxia_vlnplot
FigureS6B_Hypoxia_vlnplot
ggsave(filename = "./Figures/Epi_rename_res_045_Hallmark_HYPOXIA_module_score_vlnplot_5*3.pdf", plot = FigureS6B_Hypoxia_vlnplot, width = 5, height = 3, units = "in")

# MsigDB hallmark GLYCOLYSIS gene list
Module <- list(AddModuleScore_gene_list$MsigDB_glycolysis)
Epi_rename <- AddModuleScore(object = Epi_rename, features = Module, name = "Hallmark_GLYCOLYSIS")
FeaturePlot_scCustom(Epi_rename, features = "Hallmark_GLYCOLYSIS1", pt.size = 0.05, reduction = "umap.wnn", colors_use = c("#BBBBBB","#CCCCCC","#DDDDDD","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D"), na_cutoff = NULL) -> FigureS6B_Glycolysis
FigureS6B_Glycolysis
ggsave(filename = "./Figures/Epi_rename_res_045_Module_Hallmark_GLYCOLYSIS_wnnUMAP_4.5*4.pdf", plot = FigureS6B_Glycolysis, width = 4.5, height = 4, units = "in")

VlnPlot(Epi_rename, features = "Hallmark_GLYCOLYSIS1", pt.size = 0, cols = color.Epi.rename) + geom_boxplot(outlier.size = 0, colour = "black", fill = "white", na.rm = F, width = 0.2) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() -> FigureS6B_Glycolysis_vlnplot
FigureS6B_Glycolysis_vlnplot
ggsave(filename = "./Figures/Epi_rename_res_045_Hallmark_GLYCOLYSIS_module_score_vlnplot_5*3.pdf", plot = FigureS6B_Glycolysis_vlnplot, width = 5, height = 3, units = "in")



# Gene expression (hypoxia-induced genes)
DefaultAssay(Epi_rename) <- "RNA"

FeaturePlot_scCustom(seurat_object = Epi_rename, features = c("ADM"), pt.size = 0.05, reduction = "umap.wnn", colors_use = (brewer.pal(n=9,"OrRd"))) -> Figure6J_ADM
Figure6J_ADM
ggsave(filename = "./Figures/Epi_rename_res_045_ADM_exp_featureplot_4.5*4.pdf", plot = Figure6J_ADM, width = 4.5, height = 4, units = "in")

FeaturePlot_scCustom(seurat_object = Epi_rename, features = c("VEGFA"), pt.size = 0.05, reduction = "umap.wnn", colors_use = (brewer.pal(n=9,"OrRd"))) -> Figure6J_VEGFA
Figure6J_VEGFA
ggsave(filename = "./Figures/Epi_rename_res_045_VEGFA_exp_featureplot_4.5*4.pdf", plot = Figure6J_VEGFA, width = 4.5, height = 4, units = "in")

FeaturePlot_scCustom(seurat_object = Epi_rename, features = c("SLC2A1"), pt.size = 0.05, reduction = "umap.wnn", colors_use = (brewer.pal(n=9,"OrRd"))) -> Figure6J_SLC2A1
Figure6J_SLC2A1
ggsave(filename = "./Figures/Epi_rename_res_045_SLC2A1_exp_featureplot_4.5*4.pdf", plot = Figure6J_SLC2A1, width = 4.5, height = 4, units = "in")
