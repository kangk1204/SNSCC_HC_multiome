library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Upload RDS file
snscc_filtered <- readRDS("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/RDS_files/snscc_filtered_res_08.rds")



# Subset Epithelial cell cluster
Endo <- subset(snscc_filtered, CellType == "Endothelial cells")
DimPlot(Endo, label = TRUE, reduction = "umap.wnn")



# Gene expression data processing
DefaultAssay(Endo) <- "RNA"
Endo <- FindVariableFeatures(Endo, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
Endo <- ScaleData(Endo, features=rownames(Endo))
Endo <- RunPCA(Endo)
Endo <- RunUMAP(Endo, dims = 1:30, reduction = "harmony_pca", reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DimPlot(object = Endo, label = TRUE, reduction = "umap.rna")



# DNA accessibility data processing
DefaultAssay(Endo) <- "peaks"
Endo <- FindTopFeatures(Endo, min.cutoff = 5)
Endo <- RunTFIDF(Endo)
Endo <- RunSVD(Endo)
Endo <- RunUMAP(Endo, reduction = "harmony_lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(object = Endo, label = TRUE, reduction = "umap.atac")



# Joint UMAP visualization
# build a joint neighbor graph using both assays
Endo <- FindMultiModalNeighbors(
  object = Endo,
  reduction.list = list("harmony_pca", "harmony_lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = c("RNA.weight","peaks.weight"),
  verbose = TRUE) 

# build a joint UMAP visualization
Endo <- RunUMAP(
  object = Endo,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "umap.wnn",
  reduction.key = "wnnUMAP_",
  verbose = TRUE)



# Clustering
Endo <- FindClusters(Endo, resolution = 0.6, algorithm = 3, graph.name = "wsnn")
DimPlot(Endo, reduction = "umap.wnn", pt.size = 0.2) -> FigureS6E
FigureS6E
ggsave(filename = "./Figures/Endo_res_06_wnnUMAP_3.6*3.pdf", plot = FigureS6E, width = 3.6, height = 3, units = "in")



# scDblFiner.score boxplot
library(dittoSeq)
dittoPlot(object = Endo, plots = c("jitter","boxplot"), var = "scDblFinder.score", group.by = "seurat_clusters", boxplot.width = 0.5, boxplot.color = "black", boxplot.fill = F, boxplot.lineweight = 0.3, jitter.size = 0.1) -> FigureS6F
FigureS6F
ggsave(filename = "./Figures/Endo_res_06_scDblFinder.score_boxplot_4*3.pdf", plot = FigureS6F, width = 4, height = 3, units = "in")



# Labeling
Endo_rename <- RenameIdents(Endo, "0"="EC1","1"="EC3","2"="EC2","3"="EC4","4"="EC5","5"="Other","6"="EC6")
Endo_rename$endo_rename <- Idents(Endo_rename)



# Filter other cluster
Endo_rename_filtered <- subset(x=Endo_rename, endo_rename != "Other")

saveRDS(Endo_rename_filtered, file = "./RDS_files/Endo_rename_filtered_res_06.rds")



# Downstream analysis
Endo_rename_filtered <- readRDS("./RDS_files/Endo_rename_filtered_res_06.rds")

levels(Endo_rename_filtered) <- c("EC1","EC2","EC3","EC4","EC5","EC6")
color.Endo.filtered = c("#572D0C","#9F2E07","#CD5F09","#ECAA31","#f2df08","#D5B895")
color.sample <- c("#F7C475","#7D5008")

DimPlot(Endo_rename_filtered, label = F, reduction = "umap.wnn", cols = color.Endo.filtered, pt.size = 0.2) -> Figure6C
Figure6C
ggsave(filename = "./Figures/Endo_rename_filtered_res_06_wnnUMAP_3.8*3.pdf", plot = Figure6C, width = 3.8, height = 3, units = "in")



# Cell subtype proportion
library(dittoSeq)
dittoBarPlot(Endo_rename_filtered, var = "endo_rename", group.by = "dataset", color.panel = color.Endo.filtered, var.labels.reorder = c(1,2,3,4,5,6)) -> Figure6E
Figure6E
ggsave(filename = "./Figures/Endo_rename_filtered_res_06_endo_rename_dataset_proportion_3*4.pdf", plot = Figure6E, width = 3, height = 4, units = "in")



# Gene expression
library(scCustomize)
library(RColorBrewer)
library(ggplot2)

DefaultAssay(Endo_rename_filtered) <- "RNA"

DotPlot(Endo_rename_filtered, features = c("CXCR4","DLL4","ESM1","ENG","HSPG2","VWA1","ACKR1","SELE","SELP","CCL14","CLU","CPE","POSTN","VWF","LYVE1","PDPN","PROX1","CXCL12","DKK2","GJA4","JAG1")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) + scale_colour_gradientn(colours = c("#4575B4","#73ADD1","#AAD9E9","#FEE08F","#FDAD61","#F36D43","#D73027")) -> Figure6D
Figure6D
ggsave(filename = "./Figures/Endo_rename_filtered_dotplot2_2.5*7.pdf", plot = Figure6D, width = 7, height = 2.5, units = "in")

DotPlot(Endo_rename_filtered, features = c("ANGPT2","PLVAP","COL4A1", "COL4A2", "ITGA1", "ITGA8","HIF1A","NDRG1","ACKR1","SELP","NFKBIA","ICAM1"), cols = "RdYlBu") + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_x_discrete(limits=rev) + coord_flip() -> FigureS6I
FigureS6I
ggsave(filename = "./Figures/Endo_rename_filtered_dotplot_3.8*3.2.pdf", plot = FigureS6I, width = 3.8, height = 3.2, units = "in")

Plot_Density_Joint_Only(Endo_rename_filtered, features = c("CALCRL","RAMP2"), reduction = "umap.wnn", pt.size = 0.5, custom_palette = (brewer.pal(n=9, name = "Oranges"))) -> Figure6K_CALCRL_RAMP2
Figure6K_CALCRL_RAMP2
ggsave(filename = "./Figures/Endo_rename_filtered_CALCRL_RAMP2_joint_density_Oranges_4.4*3.pdf", plot = Figure6K_CALCRL_RAMP2, width = 4.4, height = 3, units = "in")

Plot_Density_Custom(Endo_rename_filtered, features = c("FLT1"), reduction = "umap.wnn", pt.size = 0.5, custom_palette = (brewer.pal(n=9, name = "Oranges"))) -> Figure6K_FLT1
Figure6K_FLT1
ggsave(filename = "./Figures/Endo_rename_filtered_FLT1_density_Oranges_3.9*3.pdf", plot = Figure6K_FLT1, width = 3.9, height = 3, units = "in")

Plot_Density_Custom(Endo_rename_filtered, features = c("KDR"), reduction = "umap.wnn", pt.size = 0.5, custom_palette = (brewer.pal(n=9, name = "Oranges"))) -> Figure6K_KDR
Figure6K_KDR
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Endo/Figures/Endo_rename_filtered_KDR_density_Oranges_3.9*3.pdf", plot = Figure6K_KDR, width = 3.9, height = 3, units = "in")

Plot_Density_Custom(Endo_rename_filtered, features = c("DLL4"), reduction = "umap.wnn", pt.size = 0.5, custom_palette = (brewer.pal(n=9, name = "Oranges"))) -> Figure6K_DLL4
Figure6K_DLL4
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Endo/Figures/Endo_rename_filtered_DLL4_density_Oranges_3.9*3.pdf", plot = Figure6K_DLL4, width = 3.9, height = 3, units = "in")



# AddModuleScore
DefaultAssay(Endo_rename_filtered) <- "RNA"

AddModuleScore_gene_list <- c("ACVRL1","JAG1","ANGPT1","ANGPT2","CD34","CDC42","MAPK14","TYMP","EDN1","EFNB2","EGR3","EPHB4","PTK2B","FGF2","FGFR1","VEGFD","FOXC2","FLT1","FLT4","FN1","GPLD1","NR4A1","ID1","ITGA5","ITGAV","ITGB1","KDR","LOXL2","MMP14","NOTCH1","PDGFA","PDGFRB","PGF","PIK3CA","PTGS2","PTK2","ROBO1","SHC1","SRF","TAL1","TDGF1","TEK","VAV2","VEGFA","VEGFC","FGF18","NRP1","SEMA5A","RAMP2","CIB1","ESM1","JMJD6","HEY1","ADGRA2","GREM1","SRPX2","SOX18","PARVA","RNF213","E2F8","RSPO3","OTULIN","E2F7","CCBE1","BMPER","NRARP","TNFAIP6","VCAN","SPP1","CCND2","PIK3R1","STC1","JAG2") #nature comm Dissecting the immune suppressive human prostate tumor microenvironment via integrated single-cell and spatial transcriptomic analyses

Module <- list(AddModuleScore_gene_list)
Endo_rename_filtered <- AddModuleScore(object = Endo_rename_filtered, features = Module, name = "angio")
FeaturePlot(Endo_rename_filtered, features = "angio1", pt.size = 0.01, reduction = "umap.wnn", split.by = "sample") & theme(legend.position = "right") & scale_color_gradientn(colors = rev(brewer.pal(11 ,"RdYlBu")), limits = c(-0.25, 0.6)) -> Figure6H
Figure6H
ggsave(filename = "./Figures/Endo_rename_filtered_res_06_angiogenesis_nature_comm_module_score_featureplot_splitby_sample_RdYlBu_7.5*3.pdf", plot = Figure6H, width = 7.5, height = 3, units = "in")



# DEGs
library(pheatmap)

DefaultAssay(Endo_rename_filtered) <- "RNA"

# Figure6F
Endo.all.genes <- FindMarkers(Endo_rename_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "EC1", ident.2 = c("EC2","EC3","EC4","EC5","EC6"))
write.table(Endo.all.genes, file = "./Documents/Endo_rename_filtered_res_06_FindMarkers_EC1_vs_other_DEGs_log2FC_0_pvalue_no_cutoff.txt", sep = '\t')
Endo.all.da.genes <- Endo.all.genes[Endo.all.genes$p_val_adj < 0.05, ]
write.table(Endo.all.da.genes, file = "./Documents/Endo_rename_filtered_res_06_FindMarkers_EC1_vs_other_DEGs_log2FC_058_pvalue_cutoff.txt", sep = '\t')

# FigureS6G
Endo.all.genes <- FindMarkers(Endo_rename_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 0.58, ident.1 = "EC2", ident.2 = c("EC1","EC3","EC4","EC5","EC6"))
write.table(Endo.all.genes, file = "./Documents/Endo_rename_filtered_res_06_FindMarkers_EC2_vs_other_DEGs_log2FC_0_pvalue_no_cutoff.txt", sep = '\t')
Endo.all.da.genes <- Endo.all.genes[Endo.all.genes$p_val_adj < 0.05, ]
write.table(Endo.all.da.genes, file = "./Documents/Endo_rename_filtered_res_06_FindMarkers_EC2_vs_other_DEGs_log2FC_058_pvalue_cutoff.txt", sep = '\t')
