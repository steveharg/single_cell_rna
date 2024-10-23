require("biomaRt")
require(doParallel)
library(reshape2)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
require("gplots")
require(RColorBrewer)

rm(list = ls())

source('heatmap_two_cat.r')

load('../data/processed_cts.rdata')

res_pre_cells = Cells(res_merged_so[, which(res_merged_so@meta.data$time == "pre")])
res_post_cells = Cells(res_merged_so[, which(res_merged_so@meta.data$time == "post")])
nres_pre_cells = Cells(nres_merged_so[, which(nres_merged_so@meta.data$time == "pre")])
nres_post_cells = Cells(nres_merged_so[, which(nres_merged_so@meta.data$time == "post")])

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(res_merged_so)
plot2 <- LabelPoints(plot = plot1,
                     points = top10_res_rna,
                     repel = TRUE)
pdf(
  paste(
    docs_dir,
    "plots/merged/analyses/res_rna_variable_features.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
CombinePlots(plots = list(plot1, plot2))
dev.off()

plot1 <- VariableFeaturePlot(hd_merged_so)
plot2 <- LabelPoints(plot = plot1,
                     points = top10_hd_rna,
                     repel = TRUE)
pdf(
  paste(
    docs_dir,
    "plots/merged/analyses/hd_rna_variable_features.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
CombinePlots(plots = list(plot1, plot2))
dev.off()

plot1 <- VariableFeaturePlot(nres_merged_so)
plot2 <- LabelPoints(plot = plot1,
                     points = top10_nres_rna,
                     repel = TRUE)
pdf(
  paste(
    docs_dir,
    "plots/merged/analyses/nres_rna_variable_features.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# UMAP
p1 = DimPlot(res_merged_so, reduction = "umap", group.by = "patient") + ggtitle("Responder patients")
p2 = DimPlot(
  res_merged_so,
  reduction = "umap",
  group.by = "time",
  cells = union(res_pre_cells, res_post_cells)
) + ggtitle("Responders, pre vs post")
p3 = DimPlot(res_merged_so, reduction = "umap", group.by = "stimulation") + ggtitle("Responder stimulation")
p4 = DimPlot(
  res_merged_so,
  reduction = "umap",
  group.by = "rnaClusterID",
  label = TRUE,
  do.return = TRUE
) + ggtitle("Responder clusters")
pdf(
  paste(docs_dir, "plots/merged/analyses/res_umaps.pdf", sep = ""),
  width = 10,
  height = 9
)
plot_grid(p1, p2, p3, p4)
dev.off()

pdf(
  paste(
    docs_dir,
    "plots/merged/analyses/res_umap_clusters.pdf",
    sep = ""
  ),
  width = 7.5,
  height = 6
)
p4
dev.off()

p1 = DimPlot(hd_merged_so, reduction = "umap", group.by = "stimulation") + ggtitle("HD stimulation")
p2 = DimPlot(
  hd_merged_so,
  reduction = "umap",
  group.by = "rnaClusterID",
  label = TRUE,
  do.return = TRUE
) + ggtitle("HD clusters")
pdf(
  paste(docs_dir, "plots/merged/analyses/hd_umaps.pdf", sep = ""),
  width = 16,
  height = 8
)
plot_grid(p1, p2)
dev.off()

pdf(
  paste(docs_dir, "plots/merged/analyses/hd_umap_clusters.pdf", sep = ""),
  width = 7.5,
  height = 6
)
p2
dev.off()

p1 = DimPlot(nres_merged_so, reduction = "umap", group.by = "patient") + ggtitle("Responder patients")
p2 = DimPlot(
  nres_merged_so,
  reduction = "umap",
  group.by = "time",
  cells = union(nres_pre_cells, nres_post_cells)
) + ggtitle("Non-responders, pre vs post")
p3 = DimPlot(nres_merged_so, reduction = "umap", group.by = "stimulation") + ggtitle("Non-responder stimulation")
p4 = DimPlot(
  nres_merged_so,
  reduction = "umap",
  group.by = "rnaClusterID",
  label = TRUE,
  do.return = TRUE
) + ggtitle("Non-responder clusters")
pdf(
  paste(docs_dir, "plots/merged/analyses/nres_umaps.pdf", sep = ""),
  width = 10,
  height = 9
)
plot_grid(p1, p2, p3, p4)
dev.off()

pdf(
  paste(
    docs_dir,
    "plots/merged/analyses/nres_umap_clusters.pdf",
    sep = ""
  ),
  width = 7.5,
  height = 6
)
p4
dev.off()

# Find DEGs for pre vs post from merged data
# When comparing cell populations, use
res_pre_vs_post_deg = FindMarkers(
  object = res_merged_so@assays$RNA@scale.data,
  cells.1 = res_pre_cells,
  cells.2 = res_post_cells
)

res_pre_vs_post_deg = res_pre_vs_post_deg[order(-res_pre_vs_post_deg$avg_logFC), ]

out = data.frame(gene = rownames(res_pre_vs_post_deg), res_pre_vs_post_deg[, 1:ncol(res_pre_vs_post_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "res_pre_vs_post_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

nres_pre_vs_post_deg = FindMarkers(
  object = nres_merged_so@assays$RNA@scale.data,
  cells.1 = nres_pre_cells,
  cells.2 = nres_post_cells
)

nres_pre_vs_post_deg = nres_pre_vs_post_deg[order(-nres_pre_vs_post_deg$avg_logFC), ]

out = data.frame(gene = rownames(nres_pre_vs_post_deg), nres_pre_vs_post_deg[, 1:ncol(nres_pre_vs_post_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "nres_pre_vs_post_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

save.image('../data/umap_pre_enrichment.rdata')
