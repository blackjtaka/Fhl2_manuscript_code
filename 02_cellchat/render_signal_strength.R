#!/usr/bin/env Rscript
# Render signal strength (weight-based) comparisons between SED and EX.
# Read-only — does NOT touch the original mainanalysis script.
# Outputs PDFs to ./figures_check/

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(ComplexHeatmap)
  library(grid)
})

setwd("~/Desktop/project/Reha/snRNAseq_scanpy/cellchat")
out_dir <- "figures_check"
dir.create(out_dir, showWarnings = FALSE)

cellchat.ex  <- readRDS("rds/cellchat_ex_smooth.rds")
cellchat.sed <- readRDS("rds/cellchat_sed_smooth.rds")
object.list <- list(sed = cellchat.sed, ex = cellchat.ex)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# 1. compareInteractions — total count + total signal strength (weight)
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2),
                           measure = "weight")
ggsave(file.path(out_dir, "01_compareInteractions_count_weight.pdf"),
       plot = gg1 + gg2, width = 6, height = 3.5)

# 2. netVisual_diffInteraction — circle plot of differential strength (ex vs sed)
pdf(file.path(out_dir, "02_diffInteraction_circle_count_weight.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE,
                          title.name = "Differential count (ex vs sed)")
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight",
                          title.name = "Differential signal strength (ex vs sed)")
dev.off()

# 3. netVisual_heatmap — heatmap of differential strength (ex vs sed)
ht_count  <- netVisual_heatmap(cellchat)
ht_weight <- netVisual_heatmap(cellchat, measure = "weight")
pdf(file.path(out_dir, "03_diffInteraction_heatmap_count_weight.pdf"),
    width = 11, height = 5)
draw(ht_count + ht_weight, ht_gap = unit(0.6, "cm"))
dev.off()

# 4. Per-condition circle plot: signal strength (weight) within each group
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
pdf(file.path(out_dir, "04_per_condition_signal_strength_circle.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,
                   weight.scale = TRUE, label.edge = FALSE,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Signal strength - ", names(object.list)[i]))
}
dev.off()

# 5. Numeric summary table: total weight per condition, and per-cell-type row sums
total_w_sed <- sum(cellchat.sed@net$weight)
total_w_ex  <- sum(cellchat.ex@net$weight)
total_c_sed <- sum(cellchat.sed@net$count)
total_c_ex  <- sum(cellchat.ex@net$count)

cat("\n========== Total signal summary ==========\n")
cat(sprintf("SED: total count = %d   total weight = %.4f\n",
            total_c_sed, total_w_sed))
cat(sprintf("EX : total count = %d   total weight = %.4f\n",
            total_c_ex,  total_w_ex))
cat(sprintf("Δ (ex - sed): count = %+d   weight = %+.4f   weight ratio = %.3f\n",
            total_c_ex - total_c_sed,
            total_w_ex - total_w_sed,
            total_w_ex / total_w_sed))

# Per-source weight breakdown
cell_types <- levels(cellchat.sed@idents)
df_src <- data.frame(
  cell_type = cell_types,
  out_sed   = rowSums(cellchat.sed@net$weight),
  out_ex    = rowSums(cellchat.ex@net$weight),
  in_sed    = colSums(cellchat.sed@net$weight),
  in_ex     = colSums(cellchat.ex@net$weight)
)
df_src$delta_out <- df_src$out_ex - df_src$out_sed
df_src$delta_in  <- df_src$in_ex  - df_src$in_sed
write.csv(df_src,
          file.path(out_dir, "signal_strength_by_celltype.csv"),
          row.names = FALSE)
cat("\n--- per cell type weight (sorted by Δ outgoing, ex - sed) ---\n")
print(df_src[order(-df_src$delta_out), ])

cat("\n[done] outputs in ", out_dir, "/\n", sep = "")
