library(Signac)
library(Seurat)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")

snscc_filtered$CellType_sample <- paste0(snscc_filtered$CellType,"(",snscc_filtered$sample,")")
table(snscc_filtered$CellType_sample)

snscc_filtered$CellType_dataset <- paste0(snscc_filtered$CellType,"(",snscc_filtered$dataset,")")
table(snscc_filtered$CellType_dataset)

color.snscc.filtered.CellType <- c("#d42a34", "#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")
color.sample <- c("#F7C475","#7D5008")

DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")



# Motif analysis
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

# DA motifs
library(scCustomize)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

DefaultAssay(snscc_filtered_motif) <- "chromvar"
differential.activity <- FindAllMarkers(object = snscc_filtered_motif, only.pos = T, mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0.58)
differential.activity.cut <- differential.activity[differential.activity$p_val_adj < 0.05, ]
motifLookup <- differential.activity.cut$gene
motifNames <- sapply(motifLookup, function(x) motif@motif.names[[x]])
differential.activity.cut$gene_symbol <- motifNames
differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = differential.activity.cut$gene, assays = "chromvar")

write.table(differential.activity.cut, file = "./Documents/snscc_filtered_res_08_CellType_sample_DaTF_log2FC_058.txt", sep = '\t')

differential.activity.cut %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

differential.activity.cut.average <- AverageExpression(snscc_filtered_motif, features = top10$gene, assays = "chromvar")

pheatmap::pheatmap(differential.activity.cut.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, color = colorRampPalette(rev(brewer.pal(n=11,name="PiYG")))(100), labels_row = unique(top10$gene_symbol), fontsize = 5, border_color="black") -> FigureS3F
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_sample_FindAllMarkers_avg_diff_058_DaTF_top10_5*9.pdf", plot = FigureS3F, height = 9, width = 5, units = "in")



# Sorting major cell types
DefaultAssay(snscc_filtered_motif) <- "chromvar"
snscc_filtered_motif <- ScaleData(snscc_filtered_motif, features=rownames(snscc_filtered_motif))

snscc_filtered_motif$Three_CellType <- "Immune cells"
snscc_filtered_motif$Three_CellType[snscc_filtered_motif$CellType == "Epithelial cells"] <- "Epithelial cells"
snscc_filtered_motif$Three_CellType[snscc_filtered_motif$CellType == "Fibroblasts"] <- "Fibroblasts"
snscc_filtered_motif$Three_CellType[snscc_filtered_motif$CellType == "Endothelial cells"] <- "Endothelial cells"
snscc_filtered_motif_Three_CellType <- subset(snscc_filtered_motif, Three_CellType != "Immune cells")



# Identifying mean RNA expression (CellType_sample; row expression cut-off > 1)
differential.activity.cut.RNA.MeanExp <- data.frame(AverageExpression(snscc_filtered_motif, features = differential.activity.cut$gene_symbol, assays = "RNA"))
differential.activity.cut.RNA.MeanExp$gene_symbol <- rownames(differential.activity.cut.RNA.MeanExp)

differential.activity.cut.epi.hc <- differential.activity.cut[differential.activity.cut$cluster == "Epithelial cells(HC)", ]
differential.activity.cut.fibro.hc <- differential.activity.cut[differential.activity.cut$cluster == "Fibroblasts(HC)", ]
differential.activity.cut.endo.hc <- differential.activity.cut[differential.activity.cut$cluster == "Endothelial cells(HC)", ]

differential.activity.cut.epi.tumor <- differential.activity.cut[differential.activity.cut$cluster == "Epithelial cells(Tumor)", ]
differential.activity.cut.fibro.tumor <- differential.activity.cut[differential.activity.cut$cluster == "Fibroblasts(Tumor)", ]
differential.activity.cut.endo.tumor <- differential.activity.cut[differential.activity.cut$cluster == "Endothelial cells(Tumor)", ]

differential.activity.cut.RNA.MeanExp.epi.hc <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Epithelial.cells.HC. > 1, ]
differential.activity.cut.RNA.MeanExp.fibro.hc <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Fibroblasts.HC. > 1, ]
differential.activity.cut.RNA.MeanExp.endo.hc <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Endothelial.cells.HC. > 1, ]

differential.activity.cut.RNA.MeanExp.epi.tumor <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Epithelial.cells.Tumor. > 1, ]
differential.activity.cut.RNA.MeanExp.fibro.tumor <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Fibroblasts.Tumor. > 1, ]
differential.activity.cut.RNA.MeanExp.endo.tumor <- differential.activity.cut.RNA.MeanExp[differential.activity.cut.RNA.MeanExp$RNA.Endothelial.cells.Tumor. > 1, ]

differential.activity.cut.epi.hc.RNA.MeaExp.cutoff <- differential.activity.cut.epi.hc[differential.activity.cut.epi.hc$gene_symbol %in% differential.activity.cut.RNA.MeanExp.epi.hc$gene_symbol, ]
differential.activity.cut.fibro.hc.RNA.MeaExp.cutoff <- differential.activity.cut.fibro.hc[differential.activity.cut.fibro.hc$gene_symbol %in% differential.activity.cut.RNA.MeanExp.fibro.hc$gene_symbol, ]
differential.activity.cut.endo.hc.RNA.MeaExp.cutoff <- differential.activity.cut.endo.hc[differential.activity.cut.endo.hc$gene_symbol %in% differential.activity.cut.RNA.MeanExp.endo.hc$gene_symbol, ]

differential.activity.cut.epi.tumor.RNA.MeaExp.cutoff <- differential.activity.cut.epi.tumor[differential.activity.cut.epi.tumor$gene_symbol %in% differential.activity.cut.RNA.MeanExp.epi.tumor$gene_symbol, ]
differential.activity.cut.fibro.tumor.RNA.MeaExp.cutoff <- differential.activity.cut.fibro.tumor[differential.activity.cut.fibro.tumor$gene_symbol %in% differential.activity.cut.RNA.MeanExp.fibro.tumor$gene_symbol, ]
differential.activity.cut.endo.tumor.RNA.MeaExp.cutoff <- differential.activity.cut.endo.tumor[differential.activity.cut.endo.tumor$gene_symbol %in% differential.activity.cut.RNA.MeanExp.endo.tumor$gene_symbol, ]

differential.activity.cut.RNA.MeaExp.cutoff <- rbind(differential.activity.cut.epi.hc.RNA.MeaExp.cutoff, differential.activity.cut.fibro.hc.RNA.MeaExp.cutoff, differential.activity.cut.endo.hc.RNA.MeaExp.cutoff, differential.activity.cut.epi.tumor.RNA.MeaExp.cutoff, differential.activity.cut.fibro.tumor.RNA.MeaExp.cutoff, differential.activity.cut.endo.tumor.RNA.MeaExp.cutoff)



# differential.activity.cut.RNA.MeaExp.cutoff (heatmap)
DefaultAssay(snscc_filtered_motif_Three_CellType) <- "RNA"
gene_name <- unique(differential.activity.cut.RNA.MeaExp.cutoff$gene_symbol)
motif.name <- ConvertMotifID(snscc_filtered_motif_Three_CellType, name = gene_name, assay = "peaks")

motif.average <- AverageExpression(snscc_filtered_motif_Three_CellType, features = motif.name, assays = "chromvar", slot = "scale.data")
pheatmap::pheatmap(motif.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, color = colorRampPalette(rev(brewer.pal(n=11,name="PiYG")))(15), labels_row = gene_name, fontsize = 5, border_color=alpha('black', .5)) -> figure

RNA.average <- AverageExpression(snscc_filtered_motif_Three_CellType, features = gene_name, assays = "RNA", slot = "scale.data")
pheatmap::pheatmap(RNA.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, color = colorRampPalette(rev(brewer.pal(n=9,name="RdYlBu")))(15), labels_row = gene_name, fontsize = 5, border_color=alpha('black', .5)) -> figure



# representative_motif (heatmap)
representative_gene <- c("FOSL2","FOS","JUNB","FOSL1","SREBF1","CEBPB","HIF1A","RUNX2","EGR1","NFATC2","PRDM1","ETS1","ETS2")
representative_gene.motif.name <- ConvertMotifID(snscc_filtered_motif_Three_CellType, name = representative_gene, assay = "peaks")

motif.average <- AverageExpression(snscc_filtered_motif_Three_CellType, features = representative_gene.motif.name, assays = "chromvar", slot = "scale.data")
pheatmap::pheatmap(motif.average[["chromvar"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, color = colorRampPalette(rev(brewer.pal(n=11,name="PiYG")))(15), labels_row = representative_gene, fontsize = 5, border_color=alpha('black', .5)) -> Figure3D_motif_activity
ggsave(filename = "./Figures/Epi_fibro_endo_DaTF_058_row_exp_cutoff_1_scaled_chromvar_custom_motif_heatmap_2.5*3.5.pdf", plot = Figure3D_motif_activity, height = 4, width = 2.5, units = "in")

RNA.average <- AverageExpression(snscc_filtered_motif_Three_CellType, features = representative_gene, assays = "RNA", slot = "scale.data")
pheatmap::pheatmap(RNA.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = T, color = colorRampPalette(rev(brewer.pal(n=11,name="RdYlBu")))(15), labels_row = representative_gene, fontsize = 5, border_color=alpha('black', .5)) -> Figure3D_gene_expression
ggsave(filename = "./Figures/Epi_fibro_endo_DaTF_058_row_exp_cutoff_1_scaled_chromvar_custom_motif_gene_exp_heatmap_2.5*3.5.pdf", plot = Figure3D_gene_expression, height = 4, width = 2.5, units = "in")



# Motif activity with gene expression (FeaturePlot)
# FOSL2
test_gene <- c("FOSL2")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(seurat_object = snscc_filtered_motif, features = motif.name, pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#CCCCCC","#DDDDDD","#ece6cc","#FBB4B9","#F768A1","#C51B8A","#7A0177","#49006A","#371047"), na_cutoff = -25, split.by = "sample") -> Figure3E_FOSL2_motif_activity
Figure3E_FOSL2_motif_activity
ggsave(filename = "./Figures/Motif_FOSL2_chromvar_featureplot_10*4.pdf", plot = Figure3E_FOSL2_motif_activity, height = 4, width = 10, units = "in")

DefaultAssay(snscc_filtered) <- "RNA"
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("FOSL2"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure3E_FOSL2_gene_expression
Figure3E_FOSL2_gene_expression
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_FOSL2_4*9.2.pdf", plot = Figure3E_FOSL2_gene_expression, height = 4, width = 9.2, units = "in")

# JUNB
test_gene <- c("JUNB")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(seurat_object = snscc_filtered_motif, features = motif.name, pt.size = 0.01, reduction = "umap.wnn", colors_use = c("#CCCCCC","#DDDDDD","#ece6cc","#FBB4B9","#F768A1","#C51B8A","#7A0177","#49006A","#371047"), na_cutoff = -20, split.by = "sample") -> Figure3E_JUNB_motif_activity
Figure3E_JUNB_motif_activity
ggsave(filename = "./Figures/Motif_JUNB_chromvar_featureplot_10*4.pdf", plot = Figure3E_JUNB_motif_activity, height = 4, width = 10, units = "in")

DefaultAssay(snscc_filtered) <- "RNA"
FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("JUNB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> Figure3E_JUNB_gene_expression
Figure3E_JUNB_gene_expression
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_JUNB_4*9.pdf", plot = Figure3E_JUNB_gene_expression, height = 4, width = 9, units = "in")

# NFIB
test_gene <- c("NFIB")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(seurat_object = snscc_filtered_motif, features = motif.name, pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(brewer.pal(11,"PiYG")), na_cutoff = -7.3, split.by = "sample") -> FigureS4E_NFIB_motif_activity
FigureS4E_NFIB_motif_activity
ggsave(filename = "./Figures/Motif_NFIB_chromvar_featureplot_10*4.pdf", plot = FigureS4E_NFIB_motif_activity, height = 4, width = 10, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIB"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS4E_NFIB_gene_expression
FigureS4E_NFIB_gene_expression
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFIB_4*9.2.pdf", plot = FigureS4E_NFIB_gene_expression, height = 4, width = 9.2, units = "in")

# NFIX
test_gene <- c("NFIX")
DefaultAssay(snscc_filtered_motif) <- "peaks"
motif.name <- ConvertMotifID(snscc_filtered_motif, name = test_gene)
DefaultAssay(snscc_filtered_motif) <- "chromvar"
FeaturePlot_scCustom(seurat_object = snscc_filtered_motif, features = motif.name, pt.size = 0.01, reduction = "umap.wnn", colors_use = rev(brewer.pal(11,"PiYG")), na_cutoff = -6.6, split.by = "sample") -> FigureS4E_NFIX_motif_activity
FigureS4E_NFIX_motif_activity
ggsave(filename = "./Figures/Motif_NFIX_chromvar_featureplot_10*4.pdf", plot = FigureS4E_NFIX_motif_activity, height = 4, width = 10, units = "in")

FeaturePlot_scCustom(seurat_object = snscc_filtered, features = c("NFIX"), num_columns = 2, pt.size = 0.01, reduction = "umap.wnn", colors_use = (brewer.pal(9,"OrRd")), split.by = "sample") -> FigureS4E_NFIX_gene_expression
FigureS4E_NFIX_gene_expression
ggsave(filename = "./Figures/snscc_filtered_res_08_featureplot_splitby_sample_NFIX_4*9.2.pdf", plot = FigureS4E_NFIX_gene_expression, height = 4, width = 9.2, units = "in")
