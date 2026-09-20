library(schard)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 2048 * 1024^2)
setwd("~/Desktop/project/Reha/snRNAseq_scanpy/cellchat/")

# Part 0: Data input & processing and initialization of CellChat --------

ex <- schard::h5ad2seurat('~/Desktop/project/Reha/snRNAseq_scanpy/adata/All_Ex6W.h5ad')
sed <- schard::h5ad2seurat('~/Desktop/project/Reha/snRNAseq_scanpy/adata/All_SED6W.h5ad')

ex <- NormalizeData(ex)
sed <- NormalizeData(sed)

ex$samples <- ex$sample
sed$samples <- sed$sample

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)



# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling




# Ex ----------------------------------------------------------------------

cellchat <- createCellChat(object = ex, group.by = "cell_type", assay = "RNA")

# set the used database in the object
cellchat@DB <- CellChatDB.use
# Preprocessing the expression data for cell-cell communication an --------
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)

# Compute the communication probability and infer cellular communi --------
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE) #, population.size=TRUE
cellchat <- filterCommunication(cellchat, min.cells = 10)


# Infer the cell-cell communication at a signaling pathway level ----------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate the aggregated cell-cell communication network ----------------
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,"./rds/cellchat_ex_smooth.rds")

# set up SED --------------------------------------------------------------

cellchat <- createCellChat(object = sed, group.by = "cell_type", assay = "RNA")

# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication an --------
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj=PPI.mouse)

# Compute the communication probability and infer cellular communi --------
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE) #  population.size=TRUE
cellchat <- filterCommunication(cellchat, min.cells = 10)


# Infer the cell-cell communication at a signaling pathway level ----------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate the aggregated cell-cell communication network ----------------
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,"./rds/cellchat_sed_smooth.rds")

