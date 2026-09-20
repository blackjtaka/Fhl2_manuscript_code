library(schard)
library(CellChat)
library(patchwork)

data.dir <- '~/Desktop/project/Reha/snRNAseq_scanpy/cellchat/'
setwd(data.dir)


# Load CellChat object of each dataset and merge them together ------------
cellchat.ex <- readRDS("./rds/cellchat_ex_smooth.rds")
cellchat.sed <- readRDS("./rds/cellchat_sed_smooth.rds")
object.list <- list(sed = cellchat.sed,ex = cellchat.ex)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



# Compare the total number of interactions and interaction strength -------
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# Diff interaction to CM
netVisual_diffInteraction(cellchat, weight.scale = T, targets.use = 2)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", targets.use = 2)

# heat map
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

#circle plot
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


# Compare the major sources and targets in a 2D space ---------------------

## (A) Identify cell populations with significant changes in sendin --------


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

## (B) Identify the signaling changes of specific cell populations --------
dev.new()
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cardiomyocytes")
gg1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endocardial Cells")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epicardial Cells")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3))


# Part II: Identify altered signaling with distinct network archit --------
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)



# (A) Compare the overall information flow of each signaling pathw --------
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

# (B) Compare outgoing (or incoming) signaling patterns associated --------
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



# Part III: Identify the up-gulated and down-regulated signaling l --------

#Endocardial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Epicardial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2

#CM
gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Endocardial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Epicardial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Endothelial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#FB
gg1 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Epicardial cells
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(7),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(7),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2
#Lympho EC
gg1 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Ex", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Ex", angle.x = 45, remove.isolate = T)
gg1 + gg2

