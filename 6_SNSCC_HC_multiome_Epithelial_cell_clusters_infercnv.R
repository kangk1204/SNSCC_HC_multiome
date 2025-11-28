library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Upload RDS file
Epi_rename <- readRDS("/home/iphd2/Desktop/SNSCC_HC_multiome_2024/Epi/RDS_files/Epi_rename_res_045.rds")

Idents(Epi_rename) <- "CellType_epi"
levels(Epi_rename) <- c("Basal","Ciliated","Goblet","MEC","Mucous","Serous","Tuft","TC1","TC2","TC3","TC4","TC5")

color.Epi.rename <- c("#86B6CB","#277AA5","#194F6F","#C6B1CF","#89608E","#757C98","#B3B5AA","#78010A","#E10E07","#EC796E","#E3989C","#FBC8BB")
color.sample <- c("#F7C475","#7D5008")



# infercnv
# data matrix
library(infercnv)
DefaultAssay(Epi_rename) <- "RNA"

counts_matrix <- Epi_rename@assays$RNA@counts
meta <- data.frame(labels = Idents(Epi_rename), row.names = names(Idents(Epi_rename)))
unique(meta$labels) # check the cell labels

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

gene_order <- annotation[!duplicated(annotation$gene_name),]
gene_order <- as.data.frame(gene_order[gene_order$gene_name %in% row.names(Epi_rename),])
gene_order <- gene_order[c("gene_name","seqnames","start","end")]
chrorder <- paste0("chr",c(1:22,"X","Y","MT"))
gene_order$seqnames <- factor(gene_order$seqnames, levels = chrorder)
gene_order <- with(gene_order, gene_order[order(seqnames,start),])

write.table(gene_order, file = "./infercnv/infercnv_Epi_rename_CellType_epi_gene_order.txt", sep = "\t",col.names = F, row.names = F, quote = F)



# Create the infercnv object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix, annotations_file = meta, delim = "\t", gene_order_file = "./infercnv/infercnv_Epi_rename_CellType_epi_gene_order.txt", ref_group_names = c("Ciliated","Serous","Mucous")) # analysis_mode="subclusters"
which(is.na(infercnv_obj@gene_order$chr))



# perform infercnv operations to reveal cnv signal
infercnv_obj <- infercnv::run(infercnv_obj = infercnv_obj, cutoff = 0.1, out_dir = "./infercnv/infercnv_Epi_rename_CellType_epi", denoise = T, cluster_by_groups = T, HMM = T)

saveRDS(infercnv_obj, file = "./infercnv/infercnv_Epi_rename_CellType_epi.rds")



# Plotting
infercnv_obj <- readRDS("./infercnv/infercnv_Epi_rename_CellType_epi.rds")

cluster.info <- FetchData(Epi_rename, c("CellType_epi"))

library(limma)
smoothed=apply(infercnv_obj@expr.data,2,tricubeMovingAverage,span=0.01)
cnsig=sqrt(apply((smoothed-1)^2,2,mean))
names(cnsig)=colnames(Epi_rename)

Epi_rename_cnv <- AddMetaData(Epi_rename, metadata = cnsig, col.name = "copynumber")
Epi_rename_cnv = SetIdent(Epi_rename_cnv, value = "CellType_epi")

Epi_rename$copynumber <- Epi_rename_cnv$copynumber

FeaturePlot(Epi_rename, "copynumber", pt.size = 0.01, label = F, repel = T, max.cutoff = "q95", reduction = "umap.wnn")  & scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) -> Figure5F
Figure5F
ggsave(filename = "./infercnv/infercnv_Epi_rename_copynumber_featureplot2_6.8*6.pdf", plot = Figure5F, width = 6.8, height = 6, units = "in")



# cluster reorder - inferCNV heatmap
infercnv_obj <- readRDS("./infercnv/infercnv_Epi_rename_CellType_epi.rds")

grouped_cell_indices <- infercnv_obj@observation_grouped_cell_indices
new_cluster_order <- c("TC5","TC4","TC3","TC2","TC1","Tuft","MEC","Goblet","Basal")
new_grouped_cell_indices <- grouped_cell_indices[new_cluster_order]

infercnv_obj@observation_grouped_cell_indices <- new_grouped_cell_indices

setwd("~./infercnv/infercnv_Epi_rename_CellType_epi_reorder")
infercnv::plot_cnv(infercnv_obj, output_filename = "reordered_heatmap", output_format = "png") # Figure5E
