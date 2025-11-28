library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

# Load RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")

Epi_rename <- readRDS("./RDS_files/Epi_rename_res_045.rds")

snscc_filtered$Epi <- Epi_rename$CellType_epi
snscc_filtered$CellType_Epi <- paste0(snscc_filtered$CellType,"_",snscc_filtered$Epi)
Idents(snscc_filtered) <- "CellType_Epi"
snscc_filtered <- RenameIdents(snscc_filtered, "Epithelial cells_Basal"="Basal", "Epithelial cells_Ciliated"="Ciliated", "Epithelial cells_Goblet"="Goblet", "Epithelial cells_MEC"="MEC", "Epithelial cells_Mucous"="Mucous", "Epithelial cells_Serous"="Serous", "Epithelial cells_Tuft"="Tuft", "Epithelial cells_TC1"="TC1", "Epithelial cells_TC2"="TC2", "Epithelial cells_TC3"="TC3", "Epithelial cells_TC4"="TC4", "Epithelial cells_TC5"="TC5","Fibroblasts_NA"="Fibro", "Endothelial cells_NA"="Endo", "Myeloid cells_NA"="Mye", "T cells_NA"="T", "Plasma cells_NA"="PC")
snscc_filtered$CellType_Epi <- Idents(snscc_filtered)

color.CellType_Epi <- c("#86B6CB","#277AA5","#194F6F","#C6B1CF","#89608E","#757C98","#B3B5AA","#78010A","#E10E07","#EC796E","#E3989C","#FBC8BB","#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")



# CellChat
library(CellChat)
Idents(snscc_filtered) <- "CellType_Epi"
table(snscc_filtered$CellType_Epi)
DefaultAssay(snscc_filtered) <- "RNA"
data.input <- GetAssayData(snscc_filtered, assay = "RNA", slot = "data") # normalized data matrix
meta <- snscc_filtered@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType_Epi")
summary(cellchat)
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "./cellchat/custom/CellType_epi/cellchat_snscc_filtered_res_08_CellType_Epi_secreted_signaling.rds")


# Downstream analysis
groupSize <- as.numeric(table(cellchat@idents))

# Figure6B
pdf(file = "./cellchat/custom/CellType_epi/netVisual_bubble_source_all_epi_target_endo_CALCR_MIF_VEGF_2.8*4.pdf", height = 2.8, width = 4)
netVisual_bubble(cellchat, remove.isolate = F, sources.use = c(1:12), targets.use = 14, signaling = c("CALCR","MIF","VEGF")) + scale_y_discrete(limits=rev)
dev.off()

# FigureS6D
pdf(file = "./cellchat/custom/CellType_epi/netVisual_aggregate_CALCR_5*5.pdf", height = 5, width = 5)
netVisual_aggregate(cellchat, signaling = "CALCR", layout = "circle", color.use = color.CellType_Epi)
dev.off()

pdf(file = "./cellchat/custom/CellType_epi/netVisual_aggregate_VEGF_5*5.pdf", height = 5, width = 5)
netVisual_aggregate(cellchat, signaling = "VEGF", layout = "circle", color.use = color.CellType_Epi)
dev.off()

pdf(file = "./cellchat/custom/CellType_epi/netVisual_aggregate_MIF_5*5.pdf", height = 5, width = 5)
netVisual_aggregate(cellchat, signaling = "MIF", layout = "circle", color.use = color.CellType_Epi)
dev.off()

pdf(file = "./cellchat/custom/CellType_epi/netVisual_aggregate_ncWNT_5*5.pdf", height = 5, width = 5)
netVisual_aggregate(cellchat, signaling = "ncWNT", layout = "circle", color.use = color.CellType_Epi)
dev.off()
