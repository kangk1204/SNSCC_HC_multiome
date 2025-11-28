library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)



# Load RDS file
snscc_filtered <- readRDS("./RDS_files/snscc_filtered_res_08.rds")

color.snscc.filtered.CellType <- c("#d42a34", "#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")
color.sample <- c("#F7C475","#7D5008")

snscc_filtered$CellType_sample <- paste0(snscc_filtered$CellType,"(",snscc_filtered$sample,")")
table(snscc_filtered$CellType_sample)

snscc_filtered_hc <- subset(snscc_filtered, sample == "HC")
snscc_filtered_tumor <- subset(snscc_filtered, sample == "Tumor")



# CellChatDB
library(CellChat)
options(stringsAsFactors = FALSE)

CellChatDB <- CellChatDB.human # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
write.csv(interaction_input, file = "./cellchat/custom/interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "./cellchat/custom/complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "./cellchat/custom/cofactor_input_CellChatDB.csv")
write.csv(geneInfo, file = "./cellchat/custom/geneInfo_input_CellChatDB.csv")



# Resign interaction_input_files
interaction_input <- read.csv(file = './cellchat/custom/interaction_input_CellChatDB_custom.csv', row.names = 1)
complex_input <- read.csv(file = './cellchat/custom/complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = './cellchat/custom/cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = './cellchat/custom/geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo



# CellChat
# (1) HC
library(CellChat)
Idents(snscc_filtered_hc) <- "CellType"
table(snscc_filtered_hc$CellType)
data.input <- GetAssayData(snscc_filtered_hc, assay = "RNA", slot = "data") # normalized data matrix
meta <- snscc_filtered_hc@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
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

saveRDS(cellchat, "./cellchat/custom/HC_Tumor_CellType_merge/cellchat_snscc_filtered_res_08_HC.rds")



# (2) Tumor
Idents(snscc_filtered_tumor) <- "CellType"
table(snscc_filtered_tumor$CellType)
data.input <- GetAssayData(snscc_filtered_tumor, assay = "RNA", slot = "data") # normalized data matrix
meta <- snscc_filtered_tumor@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "CellType")
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

saveRDS(cellchat, "./cellchat/custom/HC_Tumor_CellType_merge/cellchat_snscc_filtered_res_08_Tumor.rds")



# Cell-cell interaction analysis
# Merge cellchat objects
library(CellChat)
cellchat_HC <- readRDS("./cellchat/custom/HC_Tumor_CellType_merge/cellchat_snscc_filtered_res_08_HC.rds")
cellchat_Tumor <- readRDS("./cellchat/custom/HC_Tumor_CellType_merge/cellchat_snscc_filtered_res_08_Tumor.rds")
color.snscc.filtered.CellType <- c("#d42a34", "#f78c37", "#ffd827", "#62ca50", "#0677ba", "#6a32a5")
color.sample <- c("#F7C475","#7D5008")

cellchat_HC <- netAnalysis_computeCentrality(cellchat_HC, slot.name = "netP")
cellchat_Tumor <- netAnalysis_computeCentrality(cellchat_Tumor, slot.name = "netP")

object.list <- list(HC = cellchat_HC, Tumor = cellchat_Tumor)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat



# Downstream analysis
# Figure6A
pdf(file = "./cellchat/custom/HC_Tumor_CellType_merge/netVisual_circle_interaction_strength_HC_vs_Tumor_8*5.pdf", height = 5, width = 8)
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]), color.use = color.snscc.filtered.CellType, vertex.label.color = "transparent")
}
dev.off()

# FigureS6C
pdf(file = "./cellchat/custom/CellType_epi/rankNet_HC_vs_Tumor_source_Epi_target_Endo_5*3.pdf", height = 3, width = 5)
rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE, sources.use = "Epithelial cells", targets.use = "Endothelial cells", color.use = color.sample) + theme_linedraw(base_rect_size = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
