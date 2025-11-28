library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Upload RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")

snscc_filtered$CellType_sample <- paste0(snscc_filtered$CellType,"(",snscc_filtered$sample,")")
table(snscc_filtered$CellType_sample)

snscc_filtered$CellType_dataset <- paste0(snscc_filtered$CellType,"(",snscc_filtered$dataset,")")
table(snscc_filtered$CellType_dataset)

Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")



# LinkstoPeaks
DefaultAssay(snscc_filtered) <- "peaks"
snscc_filtered_regionstats <- RegionStats(snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

LinkstoPeaks <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA")
LinkstoPeaks_links <- Links(object = LinkstoPeaks)

LinkstoPeaks_links_score <- LinkstoPeaks_links[LinkstoPeaks_links$score > 0 & LinkstoPeaks_links$pvalue < 0.05, ]
write.table(LinkstoPeaks_links_score, file = "./Documents/snscc_filtered_res_08_LinkstoPeaks_links_score.txt", sep = "\t")



# DEGs (CellType_sample;FindMarkers)
DefaultAssay(snscc_filtered) <- "RNA"

snscc_filtered_Epi.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "Epithelial cells(Tumor)", ident.2 = "Epithelial cells(HC)")
snscc_filtered_Epi.all.da.genes <- snscc_filtered_Epi.all.genes[snscc_filtered_Epi.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Epi.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_Epi_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_Fibro.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "Fibroblasts(Tumor)", ident.2 = "Fibroblasts(HC)")
snscc_filtered_Fibro.all.da.genes <- snscc_filtered_Fibro.all.genes[snscc_filtered_Fibro.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Fibro.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_Endo.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "Endothelial cells(Tumor)", ident.2 = "Endothelial cells(HC)")
snscc_filtered_Endo.all.da.genes <- snscc_filtered_Endo.all.genes[snscc_filtered_Endo.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Endo.all.da.genes, file = "./Documents/snscc_filtered_res_8_CellType_sample_Endo_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_Mye.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "Myeloid cells(Tumor)", ident.2 = "Myeloid cells(HC)")
snscc_filtered_Mye.all.da.genes <- snscc_filtered_Mye.all.genes[snscc_filtered_Mye.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_Mye.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_Mye_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_T.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "T cells(Tumor)", ident.2 = "T cells(HC)")
snscc_filtered_T.all.da.genes <- snscc_filtered_T.all.genes[snscc_filtered_T.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_T.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_T_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

snscc_filtered_PC.all.genes <- FindMarkers(snscc_filtered, min.pct = 0.25, only.pos = F, logfc.threshold = 1, ident.1 = "Plasma cells(Tumor)", ident.2 = "Plasma cells(HC)")
snscc_filtered_PC.all.da.genes <- snscc_filtered_PC.all.genes[snscc_filtered_PC.all.genes$p_val_adj < 0.05, ]
write.table(snscc_filtered_PC.all.da.genes, file = "./Documents/snscc_filtered_res_08_CellType_sample_PC_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

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

snscc_filtered_res_08_FindMarkers_DEGs_log2FC_1 <- rbind(snscc_filtered_Epi.all.genes,snscc_filtered_Fibro.all.genes,snscc_filtered_Endo.all.genes,snscc_filtered_Mye.all.genes,snscc_filtered_T.all.genes,snscc_filtered_PC.all.genes)
write.table(snscc_filtered_res_08_FindMarkers_DEGs_log2FC_1, file = "./Documents/snscc_filtered_res_08_CellType_sample_FindMarkers_DEGs_log2FC_1.txt", sep = "\t")



# DARs (CellType_sample; FindMarkers)
DefaultAssay(snscc_filtered) <- "peaks"

# TSS / Distal
# split peaks between promoter regions and distal elements in each CellType
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
tss.positions <- GetTSSPositions(ranges = annotation)
tss.positions <- Extend(x=tss.positions, upstream = 1000, downstream = 100, from.midpoint = T) # obtain region around TSS

snscc_filtered_Epi.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Epithelial cells(Tumor)", ident.2 = "Epithelial cells(HC)")
snscc_filtered_Epi.all.da.peaks <- snscc_filtered_Epi.all.peaks[snscc_filtered_Epi.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Epi.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_epi <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Epi.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_epi, file = "./Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Fibro.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Fibroblasts(Tumor)", ident.2 = "Fibroblasts(HC)")
snscc_filtered_Fibro.all.da.peaks <- snscc_filtered_Fibro.all.peaks[snscc_filtered_Fibro.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Fibro.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_fibro <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Fibro.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_fibro, file = "./Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Endo.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Endothelial cells(Tumor)", ident.2 = "Endothelial cells(HC)")
snscc_filtered_Endo.all.da.peaks <- snscc_filtered_Endo.all.peaks[snscc_filtered_Endo.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Endo.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_endo <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Endo.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_endo, file = "./Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_Mye.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Myeloid cells(Tumor)", ident.2 = "Myeloid cells(HC)")
snscc_filtered_Mye.all.da.peaks <- snscc_filtered_Mye.all.peaks[snscc_filtered_Mye.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_Mye.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_mye <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_Mye.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_mye, file = "./Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_T.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "T cells(Tumor)", ident.2 = "T cells(HC)")
snscc_filtered_T.all.da.peaks <- snscc_filtered_T.all.peaks[snscc_filtered_T.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_T.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_t <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_T.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_t, file = "./Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")

snscc_filtered_PC.all.peaks <- FindMarkers(snscc_filtered, min.pct = 0.02, only.pos = F, logfc.threshold = 1, test.use = 'LR', latent.vars = 'nCount_peaks', ident.1 = "Plasma cells(Tumor)", ident.2 = "Plasma cells(HC)")
snscc_filtered_PC.all.da.peaks <- snscc_filtered_PC.all.peaks[snscc_filtered_PC.all.peaks$p_val_adj < 0.05, ]
write.table(snscc_filtered_PC.all.da.peaks, file = "./Documents/snscc_filtered_res_08_CellType_sample_PC_DARs_log2FC_1_minpct_002.txt", sep = '\t')

closest_tss_snscc_filtered_pc <- ClosestFeature(snscc_filtered, rownames(snscc_filtered_PC.all.peaks), annotation = tss.positions)
write.table(closest_tss_snscc_filtered_pc, file = "./Documents/snscc_filtered_res_08_CellType_sample_PC_DARs_log2FC_1_minpct_002_closest_tss.txt", sep = "\t")



# Making DEG files (FindMarkers)
# Epi
FindMarkers_Epi <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Epi_FindMarkers_DEGs_log2FC_058.txt", sep = '\t', header = T)

FindMarkers_Epi_N <- FindMarkers_Epi[FindMarkers_Epi$avg_log2FC < -1, ]
FindMarkers_Epi_N$cluster <- "Epi_N"
FindMarkers_Epi_T <- FindMarkers_Epi[FindMarkers_Epi$avg_log2FC > 1, ]
FindMarkers_Epi_T$cluster <- "Epi_T"
FindMarkers_Epi_N$gene <- rownames(FindMarkers_Epi_N)
FindMarkers_Epi_T$gene <- rownames(FindMarkers_Epi_T)

# Fibro
FindMarkers_Fibro <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Fibro_FindMarkers_DEGs_log2FC_1.txt", sep = '\t', header = T)

FindMarkers_Fibro_N <- FindMarkers_Fibro[FindMarkers_Fibro$avg_log2FC < -1, ]
FindMarkers_Fibro_N$cluster <- "Fibro_N"
FindMarkers_Fibro_T <- FindMarkers_Fibro[FindMarkers_Fibro$avg_log2FC > 1, ]
FindMarkers_Fibro_T$cluster <- "Fibro_T"
FindMarkers_Fibro_N$gene <- rownames(FindMarkers_Fibro_N)
FindMarkers_Fibro_T$gene <- rownames(FindMarkers_Fibro_T)

# Endo
FindMarkers_Endo <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Endo_FindMarkers_DEGs_log2FC_1.txt", sep = '\t', header = T)

FindMarkers_Endo_N <- FindMarkers_Endo[FindMarkers_Endo$avg_log2FC < -1, ]
FindMarkers_Endo_N$cluster <- "Endo_N"
FindMarkers_Endo_T <- FindMarkers_Endo[FindMarkers_Endo$avg_log2FC > 1, ]
FindMarkers_Endo_T$cluster <- "Endo_T"
FindMarkers_Endo_N$gene <- rownames(FindMarkers_Endo_N)
FindMarkers_Endo_T$gene <- rownames(FindMarkers_Endo_T)

# Mye
FindMarkers_Mye <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Mye_FindMarkers_DEGs_log2FC_1.txt", sep = '\t', header = T)

FindMarkers_Mye_N <- FindMarkers_Mye[FindMarkers_Mye$avg_log2FC < -1, ]
FindMarkers_Mye_N$cluster <- "Mye_N"
FindMarkers_Mye_T <- FindMarkers_Mye[FindMarkers_Mye$avg_log2FC > 1, ]
FindMarkers_Mye_T$cluster <- "Mye_T"
FindMarkers_Mye_N$gene <- rownames(FindMarkers_Mye_N)
FindMarkers_Mye_T$gene <- rownames(FindMarkers_Mye_T)

# T
FindMarkers_T <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_T_FindMarkers_DEGs_log2FC_1.txt", sep = '\t', header = T)

FindMarkers_T_N <- FindMarkers_T[FindMarkers_T$avg_log2FC < -1, ]
FindMarkers_T_N$cluster <- "T_N"
FindMarkers_T_T <- FindMarkers_T[FindMarkers_T$avg_log2FC > 1, ]
FindMarkers_T_T$cluster <- "T_T"
FindMarkers_T_N$gene <- rownames(FindMarkers_T_N)
FindMarkers_T_T$gene <- rownames(FindMarkers_T_T)

# PC
FindMarkers_PC <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_PC_FindMarkers_DEGs_log2FC_1.txt", sep = '\t', header = T)

FindMarkers_PC_N <- FindMarkers_PC[FindMarkers_PC$avg_log2FC < -1, ]
FindMarkers_PC_N$cluster <- "PC_N"
FindMarkers_PC_T <- FindMarkers_PC[FindMarkers_PC$avg_log2FC > 1, ]
FindMarkers_PC_T$cluster <- "PC_T"
FindMarkers_PC_N$gene <- rownames(FindMarkers_PC_N)
FindMarkers_PC_T$gene <- rownames(FindMarkers_PC_T)

# rbind
FindMarkers_DEGs <- rbind(FindMarkers_Epi_N,FindMarkers_Fibro_N,FindMarkers_Endo_N,FindMarkers_Mye_N,FindMarkers_T_N,FindMarkers_PC_N,FindMarkers_Epi_T,FindMarkers_Fibro_T,FindMarkers_Endo_T,FindMarkers_Mye_T,FindMarkers_T_T,FindMarkers_PC_T)
write.table(FindMarkers_DEGs, file = "./Documents/snscc_filtered_res_08_FindMarkers_DEGs_log2FC_1.txt", sep = "\t")



# Making DAR files (FindMarkers)
# Epi
FindMarkers_Epi <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_Epi_N <- FindMarkers_Epi[FindMarkers_Epi$avg_log2FC < -1, ]
FindMarkers_Epi_N$cluster <- "Epi_N"
FindMarkers_Epi_T <- FindMarkers_Epi[FindMarkers_Epi$avg_log2FC > 1, ]
FindMarkers_Epi_T$cluster <- "Epi_T"
FindMarkers_Epi_N$Region <- rownames(FindMarkers_Epi_N)
FindMarkers_Epi_T$Region <- rownames(FindMarkers_Epi_T)

# Fibro
FindMarkers_Fibro <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_Fibro_N <- FindMarkers_Fibro[FindMarkers_Fibro$avg_log2FC < -1, ]
FindMarkers_Fibro_N$cluster <- "Fibro_N"
FindMarkers_Fibro_T <- FindMarkers_Fibro[FindMarkers_Fibro$avg_log2FC > 1, ]
FindMarkers_Fibro_T$cluster <- "Fibro_T"
FindMarkers_Fibro_N$Region <- rownames(FindMarkers_Fibro_N)
FindMarkers_Fibro_T$Region <- rownames(FindMarkers_Fibro_T)

# Endo
FindMarkers_Endo <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_Endo_N <- FindMarkers_Endo[FindMarkers_Endo$avg_log2FC < -1, ]
FindMarkers_Endo_N$cluster <- "Endo_N"
FindMarkers_Endo_T <- FindMarkers_Endo[FindMarkers_Endo$avg_log2FC > 1, ]
FindMarkers_Endo_T$cluster <- "Endo_T"
FindMarkers_Endo_N$Region <- rownames(FindMarkers_Endo_N)
FindMarkers_Endo_T$Region <- rownames(FindMarkers_Endo_T)

# Mye
FindMarkers_Mye <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Mye_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_Mye_N <- FindMarkers_Mye[FindMarkers_Mye$avg_log2FC < -1, ]
FindMarkers_Mye_N$cluster <- "Mye_N"
FindMarkers_Mye_T <- FindMarkers_Mye[FindMarkers_Mye$avg_log2FC > 1, ]
FindMarkers_Mye_T$cluster <- "Mye_T"
FindMarkers_Mye_N$Region <- rownames(FindMarkers_Mye_N)
FindMarkers_Mye_T$Region <- rownames(FindMarkers_Mye_T)

# T
FindMarkers_T <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_T_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_T_N <- FindMarkers_T[FindMarkers_T$avg_log2FC < -1, ]
FindMarkers_T_N$cluster <- "T_N"
FindMarkers_T_T <- FindMarkers_T[FindMarkers_T$avg_log2FC > 1, ]
FindMarkers_T_T$cluster <- "T_T"
FindMarkers_T_N$Region <- rownames(FindMarkers_T_N)
FindMarkers_T_T$Region <- rownames(FindMarkers_T_T)

# PC
FindMarkers_PC <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_PC_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)

FindMarkers_PC_N <- FindMarkers_PC[FindMarkers_PC$avg_log2FC < -1, ]
FindMarkers_PC_N$cluster <- "PC_N"
FindMarkers_PC_T <- FindMarkers_PC[FindMarkers_PC$avg_log2FC > 1, ]
FindMarkers_PC_T$cluster <- "PC_T"
FindMarkers_PC_N$Region <- rownames(FindMarkers_PC_N)
FindMarkers_PC_T$Region <- rownames(FindMarkers_PC_T)

# rbind
FindMarkers_DARs <- rbind(FindMarkers_Epi_N,FindMarkers_Fibro_N,FindMarkers_Endo_N,FindMarkers_Mye_N,FindMarkers_T_N,FindMarkers_PC_N,FindMarkers_Epi_T,FindMarkers_Fibro_T,FindMarkers_Endo_T,FindMarkers_Mye_T,FindMarkers_T_T,FindMarkers_PC_T)
write.table(FindMarkers_DARs, file = "./Documents/snscc_filtered_res_08_FindMarkers_DARs_log2FC_1_minpct_002.txt", sep = "\t")



# gene linked peaks in DARs with DEGs (CellType_sample; FindMarkers)
# upload files
LinkstoPeaks_links_score <- read.table("./Documents/snscc_filtered_res_08_LinkstoPeaks_links_score.txt", sep = "\t")

FindMarkers_DARs <- read.table("./Documents/snscc_filtered_res_08_FindMarkers_DARs_log2FC_1_minpct_002.txt", sep = "\t")

FindMarkers_DEGs <- read.table("./Documents/snscc_filtered_res_08_FindMarkers_DEGs_log2FC_1.txt", sep = '\t')

# Making matrix
DARs <- FindMarkers_DARs
epi_normal_peaks <- DARs[DARs$cluster == "Epi_N", ]
epi_tumor_peaks <- DARs[DARs$cluster == "Epi_T", ]
fibro_normal_peaks <- DARs[DARs$cluster == "Fibro_N", ]
fibro_tumor_peaks <- DARs[DARs$cluster == "Fibro_T", ]
endo_normal_peaks <- DARs[DARs$cluster == "Endo_N", ]
endo_tumor_peaks <- DARs[DARs$cluster == "Endo_T", ]
mye_normal_peaks <- DARs[DARs$cluster == "Mye_N", ]
mye_tumor_peaks <- DARs[DARs$cluster == "Mye_T", ]
t_normal_peaks <- DARs[DARs$cluster == "T_N", ]
t_tumor_peaks <- DARs[DARs$cluster == "T_T", ]
pc_normal_peaks <- DARs[DARs$cluster == "PC_N", ]
pc_tumor_peaks <- DARs[DARs$cluster == "PC_T", ]

DEGs <- FindMarkers_DEGs
epi_normal_genes <- DEGs[DEGs$cluster == "Epi_N", ]
epi_tumor_genes <- DEGs[DEGs$cluster == "Epi_T", ]
fibro_normal_genes <- DEGs[DEGs$cluster == "Fibro_N", ]
fibro_tumor_genes <- DEGs[DEGs$cluster == "Fibro_T", ]
endo_normal_genes <- DEGs[DEGs$cluster == "Endo_N", ]
endo_tumor_genes <- DEGs[DEGs$cluster == "Endo_T", ]
mye_normal_genes <- DEGs[DEGs$cluster == "Mye_N", ]
mye_tumor_genes <- DEGs[DEGs$cluster == "Mye_T", ]
t_normal_genes <- DEGs[DEGs$cluster == "T_N", ]
t_tumor_genes <- DEGs[DEGs$cluster == "T_T", ]
pc_normal_genes <- DEGs[DEGs$cluster == "PC_N", ]
pc_tumor_genes <- DEGs[DEGs$cluster == "PC_T", ]

epi_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% epi_normal_peaks$Region, ]
epi_normal_peaks_genes.intersect <- epi_normal_peaks.intersect[epi_normal_peaks.intersect$gene %in% epi_normal_genes$gene, ]
epi_normal_peaks_genes.intersect$cluster <- "Epithelial cells(HC)"

fibro_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% fibro_normal_peaks$Region, ]
fibro_normal_peaks_genes.intersect <- fibro_normal_peaks.intersect[fibro_normal_peaks.intersect$gene %in% fibro_normal_genes$gene, ]
fibro_normal_peaks_genes.intersect$cluster <- "Fibroblasts(HC)"

endo_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% endo_normal_peaks$Region, ]
endo_normal_peaks_genes.intersect <- endo_normal_peaks.intersect[endo_normal_peaks.intersect$gene %in% endo_normal_genes$gene, ]
endo_normal_peaks_genes.intersect$cluster <- "Endothelial cells(HC)"

mye_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% mye_normal_peaks$Region, ]
mye_normal_peaks_genes.intersect <- mye_normal_peaks.intersect[mye_normal_peaks.intersect$gene %in% mye_normal_genes$gene, ]
mye_normal_peaks_genes.intersect$cluster <- "Myeloid cells(HC)"

t_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% t_normal_peaks$Region, ]
t_normal_peaks_genes.intersect <- t_normal_peaks.intersect[t_normal_peaks.intersect$gene %in% t_normal_genes$gene, ]
t_normal_peaks_genes.intersect$cluster <- "T cells(HC)"

pc_normal_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% pc_normal_peaks$Region, ]
pc_normal_peaks_genes.intersect <- pc_normal_peaks.intersect[pc_normal_peaks.intersect$gene %in% pc_normal_genes$gene, ]
pc_normal_peaks_genes.intersect$cluster <- "Plasma cells(HC)"

epi_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% epi_tumor_peaks$Region, ]
epi_tumor_peaks_genes.intersect <- epi_tumor_peaks.intersect[epi_tumor_peaks.intersect$gene %in% epi_tumor_genes$gene, ]
epi_tumor_peaks_genes.intersect$cluster <- "Epithelial cells(Tumor)"

fibro_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% fibro_tumor_peaks$Region, ]
fibro_tumor_peaks_genes.intersect <- fibro_tumor_peaks.intersect[fibro_tumor_peaks.intersect$gene %in% fibro_tumor_genes$gene, ]
fibro_tumor_peaks_genes.intersect$cluster <- "Fibroblasts(Tumor)"

endo_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% endo_tumor_peaks$Region, ]
endo_tumor_peaks_genes.intersect <- endo_tumor_peaks.intersect[endo_tumor_peaks.intersect$gene %in% endo_tumor_genes$gene, ]
endo_tumor_peaks_genes.intersect$cluster <- "Endothelial cells(Tumor)"

mye_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% mye_tumor_peaks$Region, ]
mye_tumor_peaks_genes.intersect <- mye_tumor_peaks.intersect[mye_tumor_peaks.intersect$gene %in% mye_tumor_genes$gene, ]
mye_tumor_peaks_genes.intersect$cluster <- "Myeloid cells(Tumor)"

t_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% t_tumor_peaks$Region, ]
t_tumor_peaks_genes.intersect <- t_tumor_peaks.intersect[t_tumor_peaks.intersect$gene %in% t_tumor_genes$gene, ]
t_tumor_peaks_genes.intersect$cluster <- "T cells(Tumor)"

pc_tumor_peaks.intersect <- LinkstoPeaks_links_score[LinkstoPeaks_links_score$peak %in% pc_tumor_peaks$Region, ]
pc_tumor_peaks_genes.intersect <- pc_tumor_peaks.intersect[pc_tumor_peaks.intersect$gene %in% pc_tumor_genes$gene, ]
pc_tumor_peaks_genes.intersect$cluster <- "Plasma cells(Tumor)"

All_DARs_gene_linked_peak_with_DEGs.intersect.tumor <- rbind(epi_tumor_peaks_genes.intersect,fibro_tumor_peaks_genes.intersect,endo_tumor_peaks_genes.intersect)
write.table(All_DARs_gene_linked_peak_with_DEGs.intersect.tumor, file = "./Documents/DEG_linked_peak_with_DEGs_intersect_minpct_002_log2FC_1_FindMarkers_Tumor_3CellType.txt", sep = "\t")

All_DARs_gene_linked_peak_with_DEGs.intersect.hc <- rbind(epi_normal_peaks_genes.intersect,fibro_normal_peaks_genes.intersect,endo_normal_peaks_genes.intersect)
write.table(All_DARs_gene_linked_peak_with_DEGs.intersect.hc, file = "./Documents/DEG_linked_peak_with_DEGs_intersect_minpct_002_log2FC_1_FindMarkers_HC_3CellType.txt", sep = "\t")

# Heatmap of DEG linked DARs
snscc_filtered$Three_CellType <- "Immune cells"
snscc_filtered$Three_CellType[snscc_filtered$CellType == "Epithelial cells"] <- "Epithelial cells"
snscc_filtered$Three_CellType[snscc_filtered$CellType == "Fibroblasts"] <- "Fibroblasts"
snscc_filtered$Three_CellType[snscc_filtered$CellType == "Endothelial cells"] <- "Endothelial cells"
snscc_filtered_Three_CellType <- subset(snscc_filtered, Three_CellType != "Immune cells")
Idents(snscc_filtered_Three_CellType)

All_DARs_gene_linked_peak_with_DEGs.intersect.hc <- read.table("./Documents/DEG_linked_peak_with_DEGs_intersect_minpct_002_log2FC_1_FindMarkers_HC_3CellType.txt", sep = "\t")

All_DARs_gene_linked_peak_with_DEGs.intersect.tumor <- read.table("./Documents/DEG_linked_peak_with_DEGs_intersect_minpct_002_log2FC_1_FindMarkers_Tumor_3CellType.txt", sep = "\t")

# Heatmap
library(scCustomize)
library(pheatmap)
library(RColorBrewer)
library(ArchR)

# DEGs linked DARs
All_DARs_gene_linked_peak_with_DEGs.intersect.hc.average <- AverageExpression(snscc_filtered_Three_CellType, assays = "peaks", features = All_DARs_gene_linked_peak_with_DEGs.intersect.hc$peak)
pheatmap::pheatmap(All_DARs_gene_linked_peak_with_DEGs.intersect.hc.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="BrBG")))(100), show_colnames = F) -> Figure3A_DARs_HC
ggsave(filename = "./Figures/FindMarkers_DEGs_linked_DARs_peaks_minpct_002_log2FC_1_heatmap_HC_3CellType_2.5*2.5.pdf", plot = Figure3A_DARs_HC, height = 2.5, width = 2.5, units = "in")

All_DARs_gene_linked_peak_with_DEGs.intersect.tumor.average <- AverageExpression(snscc_filtered_Three_CellType, assays = "peaks", features = All_DARs_gene_linked_peak_with_DEGs.intersect.tumor$peak)
pheatmap::pheatmap(All_DARs_gene_linked_peak_with_DEGs.intersect.tumor.average[["peaks"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=11,name="BrBG")))(100), show_colnames = F) -> Figure3A_DARs_Tumor
ggsave(filename = "./Figures/FindMarkers_DEGs_linked_DARs_peaks_minpct_002_log2FC_1_heatmap_Tumor_3CellType_2.5*2.5.pdf", plot = Figure3A_DARs_Tumor, height = 2.5, width = 2.5, units = "in")

# DARs linked DEGs
All_DARs_gene_linked_peak_with_DEGs.intersect.hc.average <- AverageExpression(snscc_filtered_Three_CellType, assays = "RNA", features = All_DARs_gene_linked_peak_with_DEGs.intersect.hc$gene)
pheatmap::pheatmap(All_DARs_gene_linked_peak_with_DEGs.intersect.hc.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100), show_colnames = F) -> Figure3A_DEGs_HC
ggsave(filename = "./Figures/FindMarkers_DEGs_linked_DARs_peaks_minpct_002_log2FC_1_heatmap_HC_RNA_3CellType_2.5*2.5.pdf", plot = Figure3A_DEGs_HC, height = 2.5, width = 2.5, units = "in")

All_DARs_gene_linked_peak_with_DEGs.intersect.tumor.average <- AverageExpression(snscc_filtered_Three_CellType, assays = "RNA", features = All_DARs_gene_linked_peak_with_DEGs.intersect.tumor$gene)
pheatmap::pheatmap(All_DARs_gene_linked_peak_with_DEGs.intersect.tumor.average[["RNA"]], scale = 'row', cluster_rows = F, cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(n=9,name="RdBu")))(100), show_colnames = F) -> Figure3A_DEGs_Tumor
ggsave(filename = "./Figures/FindMarkers_DEGs_linked_DARs_peaks_minpct_002_log2FC_1_heatmap_Tumor_RNA_3CellType_2.5*2.5.pdf", plot = Figure3A_DEGs_Tumor, height = 2.5, width = 2.5, units = "in")

# DARs (log2FC)
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")

FindMarkers_Epi <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Epi_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)
FindMarkers_Epi$cluster <- "Epithelial cells"
FindMarkers_Epi$gene <- rownames(FindMarkers_Epi)

FindMarkers_Fibro <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Fibro_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)
FindMarkers_Fibro$cluster <- "Fibroblasts"
FindMarkers_Fibro$gene <- rownames(FindMarkers_Fibro)

FindMarkers_Endo <- read.table("./Documents/snscc_filtered_res_08_CellType_sample_Endo_DARs_log2FC_1_minpct_002.txt", sep = '\t', header = T)
FindMarkers_Endo$cluster <- "Endothelial cells"
FindMarkers_Endo$gene <- rownames(FindMarkers_Endo)

All_DARs.intersect <- rbind(FindMarkers_Epi,FindMarkers_Fibro,FindMarkers_Endo)

Epi.DARs <- All_DARs.intersect[All_DARs.intersect$cluster == "Epithelial cells", ]
Epi.DARs <- Epi.DARs[order(Epi.DARs[, "avg_log2FC"], decreasing = T), ]
Fibro.DARs <- All_DARs.intersect[All_DARs.intersect$cluster == "Fibroblasts", ]
Fibro.DARs <- Fibro.DARs[order(Fibro.DARs[, "avg_log2FC"], decreasing = T), ]
Endo.DARs <- All_DARs.intersect[All_DARs.intersect$cluster == "Endothelial cells", ]
Endo.DARs <- Endo.DARs[order(Endo.DARs[, "avg_log2FC"], decreasing = T), ]

All_DARs.intersect.rearrange <- rbind(Epi.DARs,Fibro.DARs,Endo.DARs)

heatmap_data <- reshape(All_DARs.intersect.rearrange, idvar = "Region", timevar = "cluster", direction = "wide")

heatmap_data <- heatmap_data[, c("Region","avg_log2FC.Epithelial cells","avg_log2FC.Fibroblasts","avg_log2FC.Endothelial cells")]
rownames(heatmap_data) <- heatmap_data$Region

pheatmap::pheatmap(heatmap_data[, -1], cluster_rows = F, cluster_cols = F, show_rownames = F, fontsize_row = 5, na_col = "white", color = colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100)) -> FigureS3A
ggsave(filename = "./Figures/snscc_filtered_res_08_CellType_FindMarkers_DARs_Log2FC_1_minpct_002_3CellType_heatmap_4*6.pdf", plot = FigureS3A, width = 4, height = 6, units = "in")

snscc_filtered.all.da.peaks <- read.table("./Documents/snscc_filtered_res_08_FindMarkers_DARs_log2FC_1_minpct_002.txt", sep = "\t", header = T)
gene_counts <- table(snscc_filtered.all.da.peaks$cluster)
gene_counts # FigureS3A



# Track
DefaultAssay(snscc_filtered) <- "peaks"
Idents(snscc_filtered) <- "CellType_sample"
levels(snscc_filtered) <- c("Epithelial cells(HC)","Fibroblasts(HC)","Endothelial cells(HC)","Myeloid cells(HC)","T cells(HC)","Plasma cells(HC)","Epithelial cells(Tumor)","Fibroblasts(Tumor)","Endothelial cells(Tumor)","Myeloid cells(Tumor)","T cells(Tumor)","Plasma cells(Tumor)")

snscc_filtered_regionstats <- RegionStats(snscc_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

test_gene <- c("DSG3")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 30000, extend.downstream = -15000, idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)")) & scale_fill_manual(values = c("#fe8181","#cb2424")) & theme(axis.title.y = element_text(size = 9)) -> Figure3C_DSG3
Figure3C_DSG3
ggsave(filename = "./Figures/Track_Epi_DSG3_18*6.pdf", plot = Figure3C_DSG3, height = 6, width = 18, units = "in")

test_gene <- c("FAP")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = -60000, extend.downstream = 20000, idents = c("Fibroblasts(HC)","Fibroblasts(Tumor)")) & scale_fill_manual(values = c("#ffb38a","#ff6700")) & theme(axis.title.y = element_text(size = 9)) -> Figure3C_FAP
Figure3C_FAP
ggsave(filename = "./Figures/Track_Fibro_FAP_18*6.pdf", plot = Figure3C_FAP, height = 6, width = 18, units = "in")

test_gene <- c("ANGPT2")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = -45000, extend.downstream = 5000, idents = c("Endothelial cells(HC)","Endothelial cells(Tumor)")) & scale_fill_manual(values = c("#FADA5E","#DAA520")) & theme(axis.title.y = element_text(size = 9)) -> Figure3C_ANGPT2
Figure3C_ANGPT2
ggsave(filename = "./Figures/Track_Endo_ANGPT2_18*6.pdf", plot = Figure3C_ANGPT2, height = 6, width = 18, units = "in")

test_gene <- c("ENO1")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = -5000, extend.downstream = 8000, idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)")) & scale_fill_manual(values = c("#fe8181","#cb2424")) & theme(axis.title.y = element_text(size = 9)) -> Figure4I_ENO1
Figure4I_ENO1
ggsave(filename = "./Figures/Track_Epi_ENO1_18*6.pdf", plot = Figure4I_ENO1, height = 6, width = 18, units = "in")

test_gene <- c("KRT6A")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 500, extend.downstream = 6000, idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)")) & scale_fill_manual(values = c("#fe8181","#cb2424")) & theme(axis.title.y = element_text(size = 9)) -> Figure4I_KRT6A
Figure4I_KRT6A
ggsave(filename = "./Figures/Track_Epi_KRT6A_18*6.pdf", plot = Figure4I_KRT6A, height = 6, width = 18, units = "in")

test_gene <- c("ADM")
snscc_filtered_regionstats <- LinkPeaks(object = snscc_filtered_regionstats, peak.assay = "peaks", expression.assay = "RNA", genes.use = c(test_gene))
CoveragePlot(object = snscc_filtered_regionstats, region = test_gene, features = test_gene, peaks = T, annotation = T, links = T, extend.upstream = 5000, extend.downstream = 1500, idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)")) & scale_fill_manual(values = c("#fe8181","#cb2424")) & theme(axis.title.y = element_text(size = 9)) -> FigureS7G
FigureS7G
ggsave(filename = "./Figures/Track_Epi_ADM_15*7.pdf", plot = FigureS7G, height = 7, width = 15, units = "in")

# Vlnplot
DefaultAssay(snscc_filtered) <- "RNA"

VlnPlot(snscc_filtered, features = c("DSG3","ENO1","KRT17","NDRG1"), pt.size = 0, cols = c("#cb2424","#cb2424","#cb2424","#cb2424"), idents = c("Epithelial cells(HC)","Epithelial cells(Tumor)"), stack = T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + FontSize(y.title = 0, y.text = 10, x.title = 0, x.text = 0) -> FigureS3E_Epi
FigureS3E_Epi
ggsave(filename = "./Figures/Track_Epi_stack_vlnplot2_3*3.pdf", plot = FigureS3E_Epi, height = 3, width = 3, units = "in")

VlnPlot(snscc_filtered, features = c("FAP","FBLN2","SFRP2","TGFBI"), pt.size = 0, cols = c("#ff6700","#ff6700","#ff6700","#ff6700"), idents = c("Fibroblasts(HC)","Fibroblasts(Tumor)"), stack = T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + FontSize(y.title = 0, y.text = 10, x.title = 0, x.text = 0) -> FigureS3E_Fibro
FigureS3E_Fibro
ggsave(filename = "./Figures/Track_Fibro_stack_vlnplot2_3*3.pdf", plot = FigureS3E_Fibro, height = 3, width = 3, units = "in")

VlnPlot(snscc_filtered, features = c("ANGPT2","COL4A1","PRDM1","VMP1"), pt.size = 0, cols = c("#DAA520","#DAA520","#DAA520","#DAA520"), idents = c("Endothelial cells(HC)","Endothelial cells(Tumor)"), stack = T, flip = T) + NoLegend() + FontSize(y.title = 0, y.text = 10, x.title = 0, x.text = 0) -> FigureS3E_Endo
FigureS3E_Endo
ggsave(filename = "./Figures/Track_Endo_stack_vlnplot2_3.1*3.pdf", plot = FigureS3E_Endo, height = 3, width = 3.1, units = "in")
