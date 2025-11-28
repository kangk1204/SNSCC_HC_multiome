library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Load RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")
Epi_rename <- readRDS("./RDS_files/Epi_rename_res_045.rds")
Endo_rename_filtered <- readRDS("./RDS_files/Endo_rename_filtered_res_06.rds")

snscc_filtered$Epi <- Epi_rename$CellType_epi
snscc_filtered$Endo <- Endo_rename_filtered$endo_rename
snscc_filtered$CellType_Epi_Endo <- paste0(snscc_filtered$CellType,"_",snscc_filtered$Epi,"_",snscc_filtered$Endo)
Idents(snscc_filtered) <- "CellType_Epi_Endo"

snscc_filtered <- RenameIdents(snscc_filtered,"Endothelial cells_NA_NA"="Other")
snscc_filtered$CellType_Epi_Endo <- Idents(snscc_filtered)

snscc_filtered_Epi_Endo <- subset(x=snscc_filtered, CellType_Epi_Endo != "Other")
Idents(snscc_filtered_Epi_Endo) <- "CellType_Epi_Endo"
snscc_filtered_Epi_Endo$CellType_Epi_Endo <- Idents(snscc_filtered_Epi_Endo)
snscc_filtered_Epi_Endo <- RenameIdents(snscc_filtered_Epi_Endo, "Epithelial cells_Basal_NA"="Basal", "Epithelial cells_Ciliated_NA"="Ciliated", "Epithelial cells_Goblet_NA"="Goblet", "Epithelial cells_MEC_NA"="MEC", "Epithelial cells_Mucous_NA"="Mucous", "Epithelial cells_Serous_NA"="Serous", "Epithelial cells_Tuft_NA"="Tuft", "Epithelial cells_TC1_NA"="TC1", "Epithelial cells_TC2_NA"="TC2", "Epithelial cells_TC3_NA"="TC3", "Epithelial cells_TC4_NA"="TC4", "Epithelial cells_TC5_NA"="TC5","Fibroblasts_NA_NA"="Fibro", "Endothelial cells_NA_EC1"="EC1", "Endothelial cells_NA_EC2"="EC2", "Endothelial cells_NA_EC3"="EC3", "Endothelial cells_NA_EC4"="EC4", "Endothelial cells_NA_EC5"="EC5", "Endothelial cells_NA_EC6"="EC6", "Myeloid cells_NA_NA"="Mye", "T cells_NA_NA"="T", "Plasma cells_NA_NA"="PC")
snscc_filtered_Epi_Endo$CellType_Epi_Endo <- Idents(snscc_filtered_Epi_Endo)
Idents(snscc_filtered_Epi_Endo) <- "CellType_Epi_Endo"
levels(snscc_filtered_Epi_Endo) <- c("Basal","Ciliated","Goblet","MEC","Mucous","Serous","Tuft","TC1","TC2","TC3","TC4","TC5","Fibro","EC1","EC2","EC3","EC4","EC5","EC6","Mye","T","PC")
color.CellType_Epi_Endo_merge <- c("#86B6CB","#277AA5","#194F6F","#C6B1CF","#89608E","#757C98","#B3B5AA","#78010A","#E10E07","#EC796E","#E3989C","#FBC8BB","#f78c37","#572D0C","#9F2E07","#CD5F09","#ECAA31","#f2df08","#D5B895","#62ca50","#0677ba","#6a32a5")



# Cellchat
library(CellChat)
data.input <- GetAssayData(snscc_filtered_Epi_Endo, assay = "RNA", slot = "data") # normalized data matrix
meta <- snscc_filtered_Epi_Endo@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType_Epi_Endo")
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

saveRDS(cellchat, "./cellchat/custom/CellType_epi_endo/cellchat_snscc_filtered_res_08_Epi_Endo_res_06_merge_reorder.rds")



# Downstream analysis
groupSize <- as.numeric(table(cellchat@idents))
cellchat@netP$pathways



# CALCR signaling pathway
pathways.show <- c("CALCR")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
pairLR

# FigureS6J_ADM_CALCRL
pdf(file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/cellchat/custom/CellType_epi_endo/epi_endo_merge_netVisual_individual_ADM_CALCRL_5*5.pdf", height = 5, width = 5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[1, ], layout = "circle", color.use = color.CellType_Epi_Endo_merge)
dev.off()

# FigureS6J_ADM_RAMP2
pdf(file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/cellchat/custom/CellType_epi_endo/epi_endo_merge_netVisual_individual_ADM_RAMP2_5*5.pdf", height = 5, width = 5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[3, ], layout = "circle", color.use = color.CellType_Epi_Endo_merge)
dev.off()



# VEGF signaling pathway
pathways.show <- c("VEGF")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = F)
pairLR

# VEGFA_VEGFR1
pdf(file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/cellchat/custom/CellType_epi_endo/epi_endo_merge_netVisual_individual_VEGFA_VEGFR1_5*5.pdf", height = 5, width = 5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[1, ], layout = "circle", color.use = color.CellType_Epi_Endo_merge)
dev.off()

# FigureS6J_VEGFA_VEGFR2
pdf(file = "/home/iphd2/Desktop/SNSCC_HC_multiome_2024/cellchat/custom/CellType_epi_endo/epi_endo_merge_netVisual_individual_VEGFA_VEGFR2_5*5.pdf", height = 5, width = 5)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[2, ], layout = "circle", color.use = color.CellType_Epi_Endo_merge)
dev.off()
