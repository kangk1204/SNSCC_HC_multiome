library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# load the RNA and ATAC data
counts_HC1 <- Read10X_h5("./cellranger-arc_outputs/filtered_feature_bc_matrix.h5")
fragpath_HC1 <- "./cellranger-arc_outputs/atac_fragments.tsv.gz"

counts_HC2 <- Read10X_h5("./cellranger-arc_outputs/filtered_feature_bc_matrix.h5")
fragpath_HC2 <- "./cellranger-arc_outputs/atac_fragments.tsv.gz"

counts_Tumor1 <- Read10X_h5("./cellranger-arc_outputs/filtered_feature_bc_matrix.h5")
fragpath_Tumor1 <- "./cellranger-arc_outputs/atac_fragments.tsv.gz"

counts_Tumor2 <- Read10X_h5("./cellranger-arc_outputs/filtered_feature_bc_matrix.h5")
fragpath_Tumor2 <- "./cellranger-arc_outputs/atac_fragments.tsv.gz"

library(SoupX)
HC1_SoupX = load10X(dataDir = "./cellranger-arc_outputs/")
HC2_SoupX = load10X(dataDir = "./cellranger-arc_outputs/")
Tumor1_SoupX = load10X(dataDir = "./cellranger-arc_outputs/")
Tumor2_SoupX = load10X(dataDir = "./cellranger-arc_outputs/")

HC1_SoupX = autoEstCont(HC1_SoupX)
HC2_SoupX = autoEstCont(HC2_SoupX)
Tumor1_SoupX = autoEstCont(Tumor1_SoupX)
Tumor2_SoupX = autoEstCont(Tumor2_SoupX)

Tumor1_SoupX = setContaminationFraction(Tumor1_SoupX, 0.2)
Tumor2_SoupX = setContaminationFraction(Tumor2_SoupX, 0.2)

head(HC1_SoupX$soupProfile[order(HC1_SoupX$soupProfile$est, decreasing = T), ], n=20)
HC1.adj.matrix <- adjustCounts(HC1_SoupX)

head(HC2_SoupX$soupProfile[order(HC2_SoupX$soupProfile$est, decreasing = T), ], n=20)
HC2.adj.matrix <- adjustCounts(HC2_SoupX)

head(Tumor1_SoupX$soupProfile[order(Tumor1_SoupX$soupProfile$est, decreasing = T), ], n=20)
Tumor1.adj.matrix <- adjustCounts(Tumor1_SoupX)

head(Tumor2_SoupX$soupProfile[order(Tumor2_SoupX$soupProfile$est, decreasing = T), ], n=20)
Tumor2.adj.matrix <- adjustCounts(Tumor2_SoupX)



# create a Seurat object containing the RNA data
HC1 <- CreateSeuratObject(counts = HC1.adj.matrix, assay = "RNA")
HC2 <- CreateSeuratObject(counts = HC2.adj.matrix, assay = "RNA")
Tumor1 <- CreateSeuratObject(counts = Tumor1.adj.matrix, assay = "RNA")
Tumor2 <- CreateSeuratObject(counts = Tumor2.adj.matrix, assay = "RNA")



# Add percent.mt & percent.rb
HC1[["percent.mt"]] <- PercentageFeatureSet(HC1, pattern = "^MT-")
HC1[["percent.rb"]] <- PercentageFeatureSet(HC1, pattern = "^RP[SL]") # Ribosomal gene
HC1[["percent.hb"]] <- PercentageFeatureSet(HC1, pattern = "^HB[^(P)]") # Hemoglobin genes

HC2[["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^MT-")
HC2[["percent.rb"]] <- PercentageFeatureSet(HC2, pattern = "^RP[SL]")
HC2[["percent.hb"]] <- PercentageFeatureSet(HC2, pattern = "^HB[^(P)]")

Tumor1[["percent.mt"]] <- PercentageFeatureSet(Tumor1, pattern = "^MT-")
Tumor1[["percent.rb"]] <- PercentageFeatureSet(Tumor1, pattern = "^RP[SL]")
Tumor1[["percent.hb"]] <- PercentageFeatureSet(Tumor1, pattern = "^HB[^(P)]")

Tumor2[["percent.mt"]] <- PercentageFeatureSet(Tumor2, pattern = "^MT-")
Tumor2[["percent.rb"]] <- PercentageFeatureSet(Tumor2, pattern = "^RP[SL]")
Tumor2[["percent.hb"]] <- PercentageFeatureSet(Tumor2, pattern = "^HB[^(P)]")



# Doublet removal
library(scDblFinder)

HC1_sce = Seurat::as.SingleCellExperiment(HC1)
HC1_sce <- scDblFinder(HC1_sce)
HC1@meta.data$scDblFinder.class = HC1_sce$scDblFinder.class
HC1@meta.data$scDblFinder.score = HC1_sce$scDblFinder.score
table(HC1$scDblFinder.class)

HC2_sce = Seurat::as.SingleCellExperiment(HC2)
HC2_sce <- scDblFinder(HC2_sce)
HC2@meta.data$scDblFinder.class = HC2_sce$scDblFinder.class
HC2@meta.data$scDblFinder.score = HC2_sce$scDblFinder.score
table(HC2$scDblFinder.class)

Tumor1_sce = Seurat::as.SingleCellExperiment(Tumor1)
Tumor1_sce <- scDblFinder(Tumor1_sce)
Tumor1@meta.data$scDblFinder.class = Tumor1_sce$scDblFinder.class
Tumor1@meta.data$scDblFinder.score = Tumor1_sce$scDblFinder.score
table(Tumor1$scDblFinder.class)

Tumor2_sce = Seurat::as.SingleCellExperiment(Tumor2)
Tumor2_sce <- scDblFinder(Tumor2_sce)
Tumor2@meta.data$scDblFinder.class = Tumor2_sce$scDblFinder.class
Tumor2@meta.data$scDblFinder.score = Tumor2_sce$scDblFinder.score
table(Tumor2$scDblFinder.class)



# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))



# Creating a common peak set
# read in peak sets
peaks.HC1 <- read.table(file = "./cellranger-arc_outputs/atac_peaks.bed", col.names = c("chr", "start", "end"))
peaks.HC2 <- read.table(file = "./cellranger-arc_outputs/atac_peaks.bed", col.names = c("chr", "start", "end"))
peaks.Tumor1 <- read.table(file = "./cellranger-arc_outputs/atac_peaks.bed", col.names = c("chr", "start", "end"))
peaks.Tumor2 <- read.table(file = "./cellranger-arc_outputs/atac_peaks.bed", col.names = c("chr", "start", "end"))



# convert to genomic ranges
gr.HC1 <- makeGRangesFromDataFrame(peaks.HC1)
gr.HC2 <- makeGRangesFromDataFrame(peaks.HC2)
gr.Tumor1 <- makeGRangesFromDataFrame(peaks.Tumor1)
gr.Tumor2 <- makeGRangesFromDataFrame(peaks.Tumor2)



# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.HC1, gr.HC2, gr.Tumor1, gr.Tumor2))



# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
combined.peaks



# load metadata
md.HC1 <- read.table(
  file = "./cellranger-arc_outputs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.HC2 <- read.table(
  file = "./cellranger-arc_outputs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.Tumor1 <- read.table(
  file = "./cellranger-arc_outputs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.Tumor2 <- read.table(
  file = "./cellranger-arc_outputs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row



# perform an initial filtering of low count cells
md.HC1 <- md.HC1[md.HC1$atac_raw_reads > 1, ]
md.HC2 <- md.HC2[md.HC2$atac_raw_reads > 1, ]
md.Tumor1 <- md.Tumor1[md.Tumor1$atac_raw_reads > 1, ]
md.Tumor2 <- md.Tumor2[md.Tumor2$atac_raw_reads > 1, ]



# create fragment objects
frags.HC1 <- CreateFragmentObject(
  path = "./cellranger-arc_outputs/atac_fragments.tsv.gz",
  cells = rownames(md.HC1)
)
frags.HC2 <- CreateFragmentObject(
  path = "./cellranger-arc_outputs/atac_fragments.tsv.gz",
  cells = rownames(md.HC2)
)
frags.Tumor1 <- CreateFragmentObject(
  path = "./cellranger-arc_outputs/atac_fragments.tsv.gz",
  cells = rownames(md.Tumor1)
)
frags.Tumor2 <- CreateFragmentObject(
  path = "./cellranger-arc_outputs/atac_fragments.tsv.gz",
  cells = rownames(md.Tumor2)
)



# Quantify peaks in each dataset
HC1.counts <- FeatureMatrix(
  fragments = frags.HC1,
  features = combined.peaks,
  cells = WhichCells(HC1)
)
HC2.counts <- FeatureMatrix(
  fragments = frags.HC2,
  features = combined.peaks,
  cells = WhichCells(HC2)
)
Tumor1.counts <- FeatureMatrix(
  fragments = frags.Tumor1,
  features = combined.peaks,
  cells = WhichCells(Tumor1)
)
Tumor2.counts <- FeatureMatrix(
  fragments = frags.Tumor2,
  features = combined.peaks,
  cells = WhichCells(Tumor2)
)



# Create the objects (Add RNA assay to Seurat Object)
HC1[["ATAC"]] <- CreateChromatinAssay(counts = HC1.counts, sep = c(":", "-"), fragments = frags.HC1, annotation = annotation)
HC2[["ATAC"]] <- CreateChromatinAssay(counts = HC2.counts, sep = c(":", "-"), fragments = frags.HC2, annotation = annotation)
Tumor1[["ATAC"]] <- CreateChromatinAssay(counts = Tumor1.counts, sep = c(":", "-"), fragments = frags.Tumor1, annotation = annotation)
Tumor2[["ATAC"]] <- CreateChromatinAssay(counts = Tumor2.counts, sep = c(":", "-"), fragments = frags.Tumor2, annotation = annotation)



# Quality control
DefaultAssay(HC1) <- "ATAC"
HC1 <- NucleosomeSignal(HC1)
HC1 <- TSSEnrichment(HC1)

DefaultAssay(HC2) <- "ATAC"
HC2 <- NucleosomeSignal(HC2)
HC2 <- TSSEnrichment(HC2)

DefaultAssay(Tumor1) <- "ATAC"
Tumor1 <- NucleosomeSignal(Tumor1)
Tumor1 <- TSSEnrichment(Tumor1)

DefaultAssay(Tumor2) <- "ATAC"
Tumor2 <- NucleosomeSignal(Tumor2)
Tumor2 <- TSSEnrichment(Tumor2)

VlnPlot(object = HC1, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)
VlnPlot(object = HC2, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)
VlnPlot(object = Tumor1, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)
VlnPlot(object = Tumor2, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)



# filter out low quality cells
HC1_filtered <- subset(
  x = HC1,
  subset = nCount_ATAC < 20000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 700 &
    nCount_RNA > 1000 &
    nucleosome_signal < 0.7 &
    TSS.enrichment > 1.5 &
    percent.mt < 3 &
    scDblFinder.class == "singlet")

VlnPlot(object = snscc_filtered, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)

HC2_filtered <- subset(
  x = HC2,
  subset = nCount_ATAC < 25000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 700 &
    nCount_RNA > 1000 &
    nucleosome_signal < 0.7 &
    TSS.enrichment > 1.5 &
    percent.mt < 3 &
    scDblFinder.class == "singlet") 

VlnPlot(object = N2_filtered, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)

Tumor1_filtered <- subset(
  x = Tumor1,
  subset = nCount_ATAC < 20000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 700 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 1.5 &
    percent.mt < 3 &
    scDblFinder.class == "singlet")

VlnPlot(object = CA5_filtered, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)

Tumor2_filtered <- subset(
  x = Tumor2,
  subset = nCount_ATAC < 5000 &
    nCount_RNA < 15000 &
    nCount_ATAC > 700 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1.5 &
    percent.mt < 3 &
    scDblFinder.class == "singlet")

VlnPlot(object = CA8_filtered, features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), ncol = 5, pt.size = 0.01)



# Merge
HC1_filtered$dataset <- 'HC1'
HC2_filtered$dataset <- 'HC2'
Tumor1_filtered$dataset <- 'Tumor1'
Tumor2_filtered$dataset <- 'Tumor2'

HC1_filtered$sample <- 'HC'
HC2_filtered$sample <- 'HC'
Tumor1_filtered$sample <- 'Tumor'
Tumor2_filtered$sample <- 'Tumor'

snscc <- merge(x = HC1_filtered, y = list(HC2_filtered,Tumor1_filtered,Tumor2_filtered), add.cell.ids = c("1","2","5","8"))

snscc[["ATAC"]]

saveRDS(object = snscc, file = "./RDS_files/snscc_after_qc_before_peak_calling.rds")



#Peak calling

# call peaks using MACS2
# which macs2 in terminal - long time!!
peaks <- CallPeaks(snscc, macs2.path = "/home/anaconda3/bin/macs2", group.by = "dataset")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak - long time!!
macs2_counts <- FeatureMatrix(
  fragments = Fragments(snscc),
  features = peaks,
  cells = colnames(snscc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
snscc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(snscc),
  annotation = annotation)



# Harmony (batch correction)
# Gene expression data processing
library(harmony)
DefaultAssay(snscc) <- "RNA"
snscc <- NormalizeData(snscc, verbose = FALSE)
snscc <- FindVariableFeatures(snscc, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
snscc <- ScaleData(snscc, features=rownames(snscc))
snscc <- RunPCA(snscc)
snscc <- RunHarmony(snscc, group.by.vars = "dataset", plot_convergence = TRUE, assay.use = "RNA", reduction.save = "harmony_pca", project.dim = FALSE)
DepthCor(snscc, n = 50, reduction = "harmony_pca")
ElbowPlot(snscc, ndims = 50, reduction = "harmony_pca")
snscc <- RunUMAP(snscc, dims = 1:30, reduction = "harmony_pca", reduction.name = "umap.rna", reduction.key = "rnaUMAP_")
DimPlot(object = snscc, label = TRUE, reduction = "umap.rna")
DimPlot(object = snscc, label = TRUE, reduction = "umap.rna", group.by = "dataset")

# DNA accessibility data processing
DefaultAssay(snscc) <- "peaks"
snscc <- FindTopFeatures(snscc, min.cutoff = 5)
snscc <- RunTFIDF(snscc)
snscc <- RunSVD(snscc)
snscc <- RunHarmony(snscc, group.by.vars = "dataset", plot_convergence = TRUE, assay.use = "peaks", reduction = "lsi", reduction.save = "harmony_lsi", project.dim = FALSE)
DepthCor(snscc, n = 50, reduction = "harmony_lsi")
ElbowPlot(snscc, ndims = 50, reduction = "harmony_lsi")
snscc <- RunUMAP(snscc, reduction = "harmony_lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(object = snscc, label = TRUE, reduction = "umap.atac")
DimPlot(object = snscc, label = TRUE, reduction = "umap.atac", group.by = "dataset")

# Joint UMAP visualization
# build a joint neighbor graph using both assays
snscc <- FindMultiModalNeighbors(
  object = snscc,
  reduction.list = list("harmony_pca", "harmony_lsi"),
  modality.weight.name = c("RNA.weight","peaks.weight"),
  dims.list = list(1:30, 2:30),
  verbose = TRUE) 

# build a joint UMAP visualization
snscc <- RunUMAP(
  object = snscc,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "umap.wnn",
  reduction.key = "wnnUMAP_",
  verbose = TRUE)

saveRDS(object = snscc, file = "./RDS_files/snscc_pca_30_lsi_30_before_clustering.rds")



# Clustering
snscc <- readRDS("./RDS_files/snscc_pca_30_lsi_30_before_clustering.rds")

snscc <- FindClusters(snscc, resolution = 0.8, algorithm = 3, graph.name = "wsnn")

DimPlot(object = snscc, label = F, reduction = "umap.wnn", pt.size = 0.01) -> FigureS1D
FigureS1D
ggsave(filename = "./Figures/snscc_res_08_wnnUMAP_5*4.pdf", plot = FigureS1D, width = 5, height = 4, units = "in")



# Proportion
library(dittoSeq)
library(ggplot2)
dittoBarPlot(snscc, var = "dataset", group.by = "seurat_clusters")
dittoPlot(object = snscc, plots = c("jitter","boxplot"), var = "nFeature_RNA", group.by = "seurat_clusters", boxplot.width = 0.5, boxplot.color = "black", boxplot.fill = F, boxplot.lineweight = 0.3, jitter.size = 0.1, x.labels.rotate = F) + NoLegend() + theme(axis.text = element_text(size=17)) -> FigureS1F-a
dittoPlot(object = snscc, plots = c("jitter","boxplot"), var = "percent.mt", group.by = "seurat_clusters", boxplot.width = 0.5, boxplot.color = "black", boxplot.fill = F, boxplot.lineweight = 0.3, jitter.size = 0.1, x.labels.rotate = F) + NoLegend() + theme(axis.text = element_text(size=17)) -> FigureS1F-b
dittoPlot(object = snscc, plots = c("jitter","boxplot"), var = "scDblFinder.score", group.by = "seurat_clusters", boxplot.width = 0.5, boxplot.color = "black", boxplot.fill = F, boxplot.lineweight = 0.3, jitter.size = 0.1, x.labels.rotate = F) + NoLegend() + theme(axis.text = element_text(size=17)) -> FigureS1F-c

ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_res_08_dittoPlot_for_qc_nFeature_RNA_12*3.pdf", plot = FigureS1F-a, width = 12, height = 3, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_res_08_dittoPlot_for_qc_percent.mt_12*3.pdf", plot = FigureS1F-b, width = 12, height = 3, units = "in")
ggsave(filename = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Figures/snscc_res_08_dittoPlot_for_qc_scDblFinder.score_12*3.pdf", plot = FigureS1F-c, width = 12, height = 3, units = "in")



# Gene expression
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(scCustomize)

DefaultAssay(snscc) <- "RNA"

levels(snscc) <- c("0","1","2","3","6","7","13","14","15","16","18","19","22","24","4","5","9","8","11","23","10","12","20","17","21","25","26","27")

DotPlot(snscc, features = c("CDH1","EPCAM","KRT18","DCN","PDGFRB","THY1","CD34","CDH5","VWF","ITGAX","MS4A6A","SPI1","CD2","CD3D","CD3G","FCRL5","JCHAIN","MZB1"), cols = c("#CCCCCC","#CC0000")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+ theme_linedraw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + RotatedAxis() + scale_y_discrete(limits=rev) -> FigureS1E
FigureS1E
ggsave(filename = "./Figures/snscc_res_08_marker_gene_dotplot_6*6.pdf", plot = FigureS1E, width = 6, height = 6, units = "in")



# Rename Idents
snscc_rename <- RenameIdents(snscc, "0" = "Epithelial cells", "1" = "Epithelial cells", "2" = "Epithelial cells", "3" = "Epithelial cells", "4" = "Fibroblasts", "5" = "Fibroblasts", "6" = "Epithelial cells", "7" = "Epithelial cells", "8" = "Endothelial cells", "9" = "Fibroblasts", "10" = "Myeloid cells", "11" = "Endothelial cells", "12" = "T cells", "13" = "Epithelial cells", "14" = "Epithelial cells", "15" = "Epithelial cells", "16" = "Epithelial cells", "17" = "Other", "18" = "Epithelial cells", "19" = "Epithelial cells", "20" = "Plasma cells", "21" = "Other", "22" = "Epithelial cells", "23" = "Endothelial cells", "24" = "Epithelial cells", "25" = "Other", "26" = "Other", "27" = "Other") # res=0.8
snscc_rename$CellType <- Idents(snscc_rename)
levels(snscc_rename) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells","Other")

snscc_rename$CellType <- Idents(snscc_rename)
levels(snscc_rename) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells","Other")

DimPlot(snscc_rename, label = F, reduction = "umap.wnn", pt.size = 0.01)



# Remove low quality cells
# exept low quality cells and had no meaningful result or don't have marker gene expression
snscc_filtered <- subset(x=snscc_rename, CellType != "Other")
snscc_filtered$CellType <- Idents(snscc_filtered)
Idents(snscc_filtered) <- "CellType"
levels(snscc_filtered) <- c("Epithelial cells","Fibroblasts","Endothelial cells","Myeloid cells","T cells","Plasma cells")
color.snscc.filtered.CellType <- c("#d42a34", "#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")
DimPlot(snscc_filtered, label = TRUE, reduction = "umap.wnn", pt.size = 0.01, cols = color.snscc.filtered.CellType)
DimPlot(snscc_rename, label = TRUE, reduction = "umap.wnn", pt.size = 0.01, group.by = "dataset")

saveRDS(object = snscc_filtered, file = "./RDS_files/snscc_filtered_res_08.rds")
