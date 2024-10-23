rm(list = ls())
library(scRNAseq)
library(batchelor)
library(scater)
library(scran)
library(Seurat)
library(ggpubr)
source("shorten_gene_abseq_names.r")
source("group_enrich_individual.r")

load('../data/processed_cts.rdata')

sce_list <- list()

for (project in projects) {
  for (sample in project$samples) {
    sce <- SingleCellExperiment(assays = list(counts = sample$so@assays$RNA@counts))
    sce$sample_tag <- sample$sample_tag
    sce$stimulation <- unique(sample$so@meta.data$stimulation)
    sce$time <- unique(sample$so@meta.data$time)
    sce$patient <- unique(sample$so@meta.data$patient)
    sce$type <- unique(sample$so@meta.data$type)
    sce$project <- unique(sample$so@meta.data$project)
    altExp(sce, 'ADT') <- SingleCellExperiment(assays = list(counts = sample$so@assays$ADT@counts))
    sce_list <- c(sce_list, sce)
  }
}

p1_sce <- cbind(sce_list[[1]], sce_list[[2]], sce_list[[3]])

p2_sce <- cbind(sce_list[[4]], sce_list[[5]], sce_list[[6]], sce_list[[7]])

p3_sce <- cbind(sce_list[[8]], sce_list[[9]])

p4hd_sce <- cbind(sce_list[[10]],
                  sce_list[[11]],
                  sce_list[[12]],
                  sce_list[[13]],
                  sce_list[[14]],
                  sce_list[[15]])

p5_sce <- cbind(sce_list[[16]],
                sce_list[[17]],
                sce_list[[18]],
                sce_list[[19]],
                sce_list[[20]],
                sce_list[[21]])

p6_sce <- cbind(sce_list[[22]], sce_list[[23]])

colnames(p1_sce) <- paste0(colnames(p1_sce), "_1")
colnames(p2_sce) <- paste0(colnames(p2_sce), "_2")
colnames(p3_sce) <- paste0(colnames(p3_sce), "_3")
colnames(p4hd_sce) <- paste0(colnames(p4hd_sce), "_4")
colnames(p5_sce) <- paste0(colnames(p5_sce), "_5")
colnames(p6_sce) <- paste0(colnames(p6_sce), "_6")

p1_sce <- logNormCounts(p1_sce, use_altexps = TRUE)
p2_sce <- logNormCounts(p2_sce, use_altexps = TRUE)
p3_sce <- logNormCounts(p3_sce, use_altexps = TRUE)
p4hd_sce <- logNormCounts(p4hd_sce, use_altexps = TRUE)
p5_sce <- logNormCounts(p5_sce, use_altexps = TRUE)
p6_sce <- logNormCounts(p6_sce, use_altexps = TRUE)

dec1 <- modelGeneVar(p1_sce)
dec1.adt <- modelGeneVar(altExp(p1_sce, 'ADT'))
dec2 <- modelGeneVar(p2_sce)
dec2.adt <- modelGeneVar(altExp(p2_sce, 'ADT'))
dec3 <- modelGeneVar(p3_sce)
dec3.adt <- modelGeneVar(altExp(p3_sce, 'ADT'))
dec4 <- modelGeneVar(p4hd_sce)
# dec4.adt <- modelGeneVar(altExp(p4hd_sce, 'ADT'))
dec5 <- modelGeneVar(p5_sce)
dec5.adt <- modelGeneVar(altExp(p5_sce, 'ADT'))
dec6 <- modelGeneVar(p6_sce)
dec6.adt <- modelGeneVar(altExp(p6_sce, 'ADT'))

# Batch-correct RNA
combined.dec <- combineVar(dec1, dec2, dec3, dec4, dec5, dec6)
chosen.hvgs <- combined.dec$bio > 0
summary(chosen.hvgs)

combined <- correctExperiments(
  "Responder Project 1 (Patient 19&24)" = p1_sce,
  "Responder Project 2 (Patient 19&24 O/N)" =
    p2_sce,
  "Responder Project 3 (Patient 18)" = p3_sce,
  "Healthy donor Project 1" = p4hd_sce,
  "Non-responder Project 2 (Patient 20, 21, 22)" =
    p5_sce,
  "Non-responder Project 3 (Patient 22)" =
    p6_sce,
  PARAM = NoCorrectParam()
)

combined <- runPCA(combined, subset_row = chosen.hvgs)
pdf(
  file = paste(
    docs_dir,
    "plots/b_corrected/1.pca/1.pca_rna_no_batch_correct.pdf",
    sep = ""
  ),
  width = 8,
  height = 5
)
plotPCA(combined, colour_by = "batch")
dev.off()

f.out <- fastMNN(
  "Responder Project 1 (Patient 19&24)" = p1_sce,
  "Responder Project 2 (Patient 19&24 O/N)" = p2_sce,
  "Responder Project 3 (Patient 18)" = p3_sce,
  "Healthy donor Project 1" = p4hd_sce,
  "Non-responder Project 2 (Patient 20, 21, 22)" = p5_sce,
  "Non-responder Project 3 (Patient 22)" = p6_sce,
  subset.row = chosen.hvgs
)
str(reducedDim(f.out, "corrected"))
rle(f.out$batch)

all_colData <- rbind(
  colData(p1_sce),
  colData(p2_sce),
  colData(p3_sce),
  colData(p4hd_sce),
  colData(p5_sce),
  colData(p6_sce)
)

colData(f.out) <- cbind(colData(f.out), all_colData)

pdf(
  file = paste(
    docs_dir,
    "plots/b_corrected/1.pca/2.pca_rna_batch_corrected.pdf",
    sep = ""
  ),
  width = 8,
  height = 5
)
plotReducedDim(f.out, colour_by = "batch", dimred = "corrected")
dev.off()

# Batch-correct ADT
# Don't choose a subset of proteins, we only have 9 anyway
# combined.adt.dec <- combineVar(dec1.adt, dec2.adt, dec3.adt,
#                                dec5.adt, dec6.adt)
# chosen.adt.hvgs <- combined.adt.dec$bio > 0
# summary(chosen.adt.hvgs)

combined.adt <- correctExperiments(
  "Responder Project 1 (Patient 19&24)" = altExp(p1_sce),
  "Responder Project 2 (Patient 19&24 O/N)" =
    altExp(p2_sce),
  "Responder Project 3 (Patient 18)" =
    altExp(p3_sce),
  "Healthy donor Project 1" = altExp(p4hd_sce),
  "Non-responder Project 2 (Patient 20, 21, 22)" =
    altExp(p5_sce),
  "Non-responder Project 3 (Patient 22)" =
    altExp(p6_sce),
  PARAM = NoCorrectParam()
)
combined.adt <- runPCA(combined.adt)

pdf(
  file = paste(
    docs_dir,
    "plots/b_corrected/1.pca/3.pca_adt_no_batch_correct.pdf",
    sep = ""
  ),
  width = 8,
  height = 5
)
plotPCA(combined.adt, colour_by = "batch")
dev.off()

f.adt.out <- fastMNN(
  "Responder Project 1 (Patient 19&24)" = altExp(p1_sce),
  "Responder Project 2 (Patient 19&24 O/N)" = altExp(p2_sce),
  "Responder Project 3 (Patient 18)" = altExp(p3_sce),
  "Healthy donor Project 1" = altExp(p4hd_sce),
  "Non-responder Project 2 (Patient 20, 21, 22)" = altExp(p5_sce),
  "Non-responder Project 3 (Patient 22)" = altExp(p6_sce)
)
str(reducedDim(f.adt.out, "corrected"))
rle(f.adt.out$batch)

all_colData.adt <- rbind(
  colData(altExp(p1_sce)),
  colData(altExp(p2_sce)),
  colData(altExp(p3_sce)),
  colData(altExp(p4hd_sce)),
  colData(altExp(p5_sce)),
  colData(altExp(p6_sce))
)

colData(f.adt.out) <- cbind(colData(f.adt.out), all_colData.adt)

pdf(
  file = paste(
    docs_dir,
    "plots/b_corrected/1.pca/4.pca_adt_batch_corrected.pdf",
    sep = ""
  ),
  width = 8,
  height = 5
)
plotReducedDim(f.adt.out, colour_by = "batch", dimred = "corrected")
dev.off()

assay(f.out) <- as.matrix(assay(f.out))
b_corrected_so <- as.Seurat(f.out, data = "reconstructed", counts = NULL)

################################################################################
# Now need to add the ADT resuts as an extra assay to the seurat object,
# and see if we can plot ADT expression on the RNA UMAP
################################################################################

assay(f.adt.out) <- as.matrix(assay(f.adt.out))
b_corrected_so.adt <- as.Seurat(f.adt.out,
                                data = "reconstructed",
                                counts = NULL,
                                assay = "ADT")
b_corrected_so[["ADT"]] <- CreateAssayObject(data = assay(f.adt.out))

# RNA
b_corrected_so <- FindVariableFeatures(b_corrected_so)
b_corrected_so <- ScaleData(b_corrected_so)
b_corrected_so <- RunPCA(b_corrected_so, verbose = FALSE)

pdf(
  file = paste(
    docs_dir,
    "plots/b_corrected/1.pca/5.pca_elbow_plot.pdf",
    sep = ""
  ),
  width = 6,
  height = 5
)
ElbowPlot(b_corrected_so)
dev.off()

num_pcs <- 15

b_corrected_so <- FindNeighbors(b_corrected_so, dims = 1:num_pcs)
b_corrected_so <- FindClusters(b_corrected_so, resolution = rna_clusters_resolution)
b_corrected_so <- RunUMAP(b_corrected_so,
                          reduction = "pca",
                          dims = 1:num_pcs)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/1.umap_by_batch.pdf",
    sep = ""
  ),
  width = 14,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap", group.by = "batch") + ggtitle("UMAP by batch")
dev.off()

# ADT
DefaultAssay(b_corrected_so) <- "ADT"
b_corrected_so <- FindVariableFeatures(b_corrected_so)
b_corrected_so <- ScaleData(b_corrected_so)

# Back to default RNA
DefaultAssay(b_corrected_so) <- "RNA"

# Test

corrected_projects <- SplitObject(b_corrected_so, split.by = "project")
# Remove healthy donor project as there is no time metadata
corrected_projects <- corrected_projects[-4]

for (corrected_project in corrected_projects) {
  print(unique(corrected_project$project))
  
  file_proj_name <- gsub(
    x = unique(corrected_project$project),
    pattern = '/',
    replacement = ''
  )
  
  corrected_project_times <- SplitObject(corrected_project, split.by = "time")
  
  result = tryCatch({
    pre_vs_post_rna_b_corrected_deg = FindMarkers(
      object = b_corrected_so@assays$RNA@scale.data,
      cells.1 = Cells(corrected_project_times$pre),
      cells.2 = Cells(corrected_project_times$post),
      logfc.threshold = 0.01
    )
    
    pre_vs_post_adt_b_corrected_deg = FindMarkers(
      object = b_corrected_so@assays$ADT@scale.data,
      cells.1 = Cells(corrected_project_times$pre),
      cells.2 = Cells(corrected_project_times$post),
      logfc.threshold = 0.01
    )
    
    out = data.frame(gene = rownames(pre_vs_post_rna_b_corrected_deg),
                     pre_vs_post_rna_b_corrected_deg[, 1:ncol(pre_vs_post_rna_b_corrected_deg)])
    write.table(
      x = out,
      file = paste0(
        data_dir,
        file_proj_name,
        "_pre_vs_post_rna_b_corrected_deg.tsv"
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
    
    out = data.frame(gene = rownames(pre_vs_post_adt_b_corrected_deg),
                     pre_vs_post_adt_b_corrected_deg[, 1:ncol(pre_vs_post_adt_b_corrected_deg)])
    write.table(
      x = out,
      file = paste0(
        data_dir,
        file_proj_name,
        "_pre_vs_post_adt_b_corrected_deg.tsv"
      ),
      sep = "\t",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE
    )
    
  }, warning = function(warning_condition) {
    print(warning_condition)
  }, error = function(error_condition) {
    print(error_condition)
  }, finally = {
    print("tryCatch finally")
  })
}

# Notes on the results from the tests above:
# For RNA, results here are much more plausible than those from Seurat's method (fold-change is up/down in line with the non-batch-corrected results)
# For ADT, we get much better p-vals than we did we either the (Seurat) merged or integrated methods

# Analyse

# Extra plot not requested - UMAP with all cells, coloured by cluster
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/2.umap_all_cells_snn_clusters.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap") + ggtitle("UMAP all cells, SNN clusters")
dev.off()

# Extra plot not requested - UMAP with R, NR, HD all superimposed
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/3.umap_all_cells_by_type.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap", group.by = "type") + ggtitle("UMAP all cells, grouped by type")
dev.off()

# UMAP responders only, all-cell cluster numbers
resp_cells = Cells(b_corrected_so[, which(b_corrected_so@meta.data$type == "resp")])
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/4.umap_res_snn_clusters.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap", cells = resp_cells) + ggtitle("UMAP responders, SNN clusters")
dev.off()

# # UMAP responders only, responders cluster numbers
resp_so <- subset(x = b_corrected_so, cells = resp_cells)
# resp_so <- FindNeighbors(resp_so, dims = 1:num_pcs)
# resp_so <- FindClusters(resp_so, resolution = rna_clusters_resolution)
# resp_so <- RunUMAP(resp_so, reduction = "pca", dims = 1:num_pcs)
# DimPlot(resp_so, reduction = "umap") + ggtitle("UMAP responders, coloured by responder clusters")

# UMAP non-responders only
nresp_cells = Cells(b_corrected_so[, which(b_corrected_so@meta.data$type == "nonresp")])
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/5.nres_umap_snn_clusters.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap", cells = nresp_cells) + ggtitle("UMAP non-responders, SNN clusters")
dev.off()

# # UMAP non-responders only, non-responders cluster numbers
nresp_so <- subset(x = b_corrected_so, cells = nresp_cells)
# nresp_so <- FindNeighbors(nresp_so, dims = 1:num_pcs)
# nresp_so <- FindClusters(nresp_so, resolution = rna_clusters_resolution)
# nresp_so <- RunUMAP(nresp_so, reduction = "pca", dims = 1:num_pcs)
# DimPlot(nresp_so, reduction = "umap") + ggtitle("UMAP non-responders, coloured by non-responder clusters")

# UMAP hd only
hd_cells = Cells(b_corrected_so[, which(b_corrected_so@meta.data$type == "hd")])
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/6.hd_umap_snn_clusters.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(b_corrected_so, reduction = "umap", cells = hd_cells) + ggtitle("UMAP healthy donors, SNN clusters")
dev.off()

# # UMAP hd only, hd cluster numbers
hd_so <- subset(x = b_corrected_so, cells = hd_cells)
# hd_so <- FindNeighbors(hd_so, dims = 1:num_pcs)
# hd_so <- FindClusters(hd_so, resolution = rna_clusters_resolution)
# hd_so <- RunUMAP(hd_so, reduction = "pca", dims = 1:num_pcs)
# DimPlot(hd_so, reduction = "umap") + ggtitle("UMAP healthy donors, coloured by healthy donor clusters")

# UMAP pre vs post (regardless of whether resp or nresp)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/7.pre_vs_post_umap.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(
  b_corrected_so,
  reduction = "umap",
  cells = c(nresp_cells, resp_cells),
  group.by = "time"
) + ggtitle("UMAP pre vs post, responders and non-responders")
dev.off()

# UMAP pre vs post, responders only
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/8.pre_vs_post_resp_umap.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(
  b_corrected_so,
  reduction = "umap",
  cells = resp_cells,
  group.by = "time"
) + ggtitle("UMAP pre vs post, responders")
dev.off()

# UMAP pre vs post, non-responders only
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/2.umaps/9.pre_vs_post_nresp_umap.pdf",
    sep = ""
  ),
  width = 11,
  height = 10
)
DimPlot(
  b_corrected_so,
  reduction = "umap",
  cells = nresp_cells,
  group.by = "time"
) + ggtitle("UMAP pre vs post, non-responders")
dev.off()


# Find DEGs for pre vs post, responders
res_pre_cells = Cells(resp_so[, which(resp_so@meta.data$time == "pre")])
res_post_cells = Cells(resp_so[, which(resp_so@meta.data$time == "post")])
res_pre_vs_post_deg = FindMarkers(
  object = b_corrected_so@assays$RNA@scale.data,
  cells.1 = res_pre_cells,
  cells.2 = res_post_cells
)

out = data.frame(gene = rownames(res_pre_vs_post_deg), res_pre_vs_post_deg[, 1:ncol(res_pre_vs_post_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/1.res_pre_vs_post_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

res_pre_vs_post_adt_deg = FindMarkers(
  object = b_corrected_so@assays$ADT@scale.data,
  cells.1 = res_pre_cells,
  cells.2 = res_post_cells
)

out = data.frame(gene = rownames(res_pre_vs_post_adt_deg), res_pre_vs_post_adt_deg[, 1:ncol(res_pre_vs_post_adt_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/2.res_pre_vs_post_adt_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Find DEGs for pre vs post, non-responders
nres_pre_cells = Cells(nresp_so[, which(nresp_so@meta.data$time == "pre")])
nres_post_cells = Cells(nresp_so[, which(nresp_so@meta.data$time == "post")])
nres_pre_vs_post_deg = FindMarkers(
  object = b_corrected_so@assays$RNA@scale.data,
  cells.1 = nres_pre_cells,
  cells.2 = nres_post_cells
)

out = data.frame(gene = rownames(nres_pre_vs_post_deg), nres_pre_vs_post_deg[, 1:ncol(nres_pre_vs_post_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/3.nres_pre_vs_post_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

nres_pre_vs_post_adt_deg = FindMarkers(
  object = b_corrected_so@assays$ADT@scale.data,
  cells.1 = nres_pre_cells,
  cells.2 = nres_post_cells
)

out = data.frame(gene = rownames(nres_pre_vs_post_adt_deg), nres_pre_vs_post_adt_deg[, 1:ncol(nres_pre_vs_post_adt_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/4.nres_pre_vs_post_adt_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Find DEGs for HD vs post responders
post_cells = Cells(b_corrected_so[, which(b_corrected_so@meta.data$time == "post")])
post_resp_cells = intersect(resp_cells, post_cells)
post_nresp_cells = intersect(nresp_cells, post_cells)
hd_post_resp_so <- subset(x = b_corrected_so, cells = union(hd_cells, post_resp_cells)) ##### up to here
hd_post_resp_so$type[hd_post_resp_so$type == "resp"] <- "post resp"
hd_post_nresp_so <- subset(x = b_corrected_so, cells = union(hd_cells, post_nresp_cells))
hd_post_nresp_so$type[hd_post_nresp_so$type == "nonresp"] <- "post nonresp"

hd_vs_post_res_deg = FindMarkers(
  object = hd_post_resp_so@assays$RNA@scale.data,
  cells.1 = hd_cells,
  cells.2 = res_post_cells
)

out = data.frame(gene = rownames(hd_vs_post_res_deg), hd_vs_post_res_deg[, 1:ncol(hd_vs_post_res_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/5.hd_vs_post_res_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

hd_vs_post_res_adt_deg = FindMarkers(
  object = hd_post_resp_so@assays$ADT@scale.data,
  cells.1 = hd_cells,
  cells.2 = res_post_cells
)

out = data.frame(gene = rownames(hd_vs_post_res_adt_deg), hd_vs_post_res_adt_deg[, 1:ncol(hd_vs_post_res_adt_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/6.hd_vs_post_res_adt_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# Find DEGs for HD vs post non-responders
hd_vs_post_nres_deg = FindMarkers(
  object = hd_post_nresp_so@assays$RNA@scale.data,
  cells.1 = hd_cells,
  cells.2 = nres_post_cells
)

out = data.frame(gene = rownames(hd_vs_post_nres_deg), hd_vs_post_nres_deg[, 1:ncol(hd_vs_post_nres_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/7.hd_vs_post_nres_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

hd_vs_post_nres_adt_deg = FindMarkers(
  object = hd_post_nresp_so@assays$ADT@scale.data,
  cells.1 = hd_cells,
  cells.2 = nres_post_cells
)

out = data.frame(gene = rownames(hd_vs_post_nres_adt_deg), hd_vs_post_nres_adt_deg[, 1:ncol(hd_vs_post_nres_adt_deg)])
write.table(
  x = out,
  file = paste0(data_dir, "b_corrected/8.hd_vs_post_nres_adt_deg.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# GO terms and KEGG Pathways
enrich(res_pre_vs_post_deg,
       "1.res_pre_vs_post_deg_b_corrected",
       "Responders, pre vs post")
enrich(
  nres_pre_vs_post_deg,
  "2.nres_pre_vs_post_deg_b_corrected",
  "Non-responders, pre vs post"
)
enrich(hd_vs_post_res_deg,
       "3.hd_vs_post_res_deg_b_corrected",
       "HD vs post responders")
enrich(
  hd_vs_post_nres_deg,
  "4.hd_vs_post_nres_deg_b_corrected",
  "HD vs post non-responders"
)

# Heatmaps
nahids_genes <- c(
  "CD4",
  "CTLA4",
  "FOXP3",
  "CD69",
  "IFNG",
  "IL2",
  "CD8A",
  "HAVCR2",
  "PDCD1",
  "FAS",
  "TNFRSF9",
  "TNF",
  "IL17A",
  "IL2RA",
  "IL7R",
  "CD3D",
  "LAG3",
  "CD40LG",
  "HLA.DRA"
)

res_pre_vs_post_deg <- res_pre_vs_post_deg[res_pre_vs_post_deg$p_val_adj < 0.05, ]
res_pre_vs_post_deg <- res_pre_vs_post_deg[order(res_pre_vs_post_deg$avg_logFC), ]
rownames(res_pre_vs_post_deg)[rownames(res_pre_vs_post_deg) %in% nahids_genes]

resp_so$time <- factor(resp_so$time, levels = c("pre", "post"))
hist(resp_so[rownames(res_pre_vs_post_deg)[rownames(res_pre_vs_post_deg) %in% nahids_genes], ]@assays$RNA@scale.data)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/1.rna_heatmap_resp_pre_vs_post.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  resp_so,
  features = rownames(res_pre_vs_post_deg)[rownames(res_pre_vs_post_deg) %in% nahids_genes],
  group.by = 'time',
  size = 2.5,
  disp.min = -2,
  disp.max = 5
) +
  scale_fill_gradientn(
    limits = c(-3, 5),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Responders, RNA")
dev.off()

res_pre_vs_post_adt_deg <- res_pre_vs_post_adt_deg[res_pre_vs_post_adt_deg$p_val_adj < 0.05, ]
res_pre_vs_post_adt_deg <- res_pre_vs_post_adt_deg[order(res_pre_vs_post_adt_deg$avg_logFC), ]

hist(resp_so@assays$ADT@scale.data)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/2.adt_heatmap_resp_pre_vs_post.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  resp_so,
  features = rownames(res_pre_vs_post_adt_deg),
  group.by = 'time',
  size = 2.5,
  disp.min = -2,
  disp.max = 2,
  assay = "ADT"
) +
  scale_fill_gradientn(
    limits = c(-2, 2),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Responders, proteins")
dev.off()

nres_pre_vs_post_deg <- nres_pre_vs_post_deg[nres_pre_vs_post_deg$p_val_adj < 0.05, ]
nres_pre_vs_post_deg <- nres_pre_vs_post_deg[order(nres_pre_vs_post_deg$avg_logFC), ]
rownames(nres_pre_vs_post_deg)[rownames(nres_pre_vs_post_deg) %in% nahids_genes]

nresp_so$time <- factor(nresp_so$time, levels = c("pre", "post"))
hist(nresp_so[rownames(nres_pre_vs_post_deg)[rownames(nres_pre_vs_post_deg) %in% nahids_genes], ]@assays$RNA@scale.data)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/3.rna_heatmap_nresp_pre_vs_post.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  nresp_so,
  features = rownames(nres_pre_vs_post_deg)[rownames(nres_pre_vs_post_deg) %in% nahids_genes],
  group.by = 'time',
  size = 2.5,
  disp.min = -2,
  disp.max = 3
) +
  scale_fill_gradientn(
    limits = c(-2, 3),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Non-responders, RNA")
dev.off()

nres_pre_vs_post_adt_deg <- nres_pre_vs_post_adt_deg[nres_pre_vs_post_adt_deg$p_val_adj < 0.05, ]
nres_pre_vs_post_adt_deg <- nres_pre_vs_post_adt_deg[order(nres_pre_vs_post_adt_deg$avg_logFC), ]

hist(nresp_so@assays$ADT@scale.data[rownames(nres_pre_vs_post_adt_deg)[1:3], ])
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/4.adt_heatmap_nresp_pre_vs_post.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  nresp_so,
  features = rownames(nres_pre_vs_post_adt_deg)[1:3],
  group.by = 'time',
  size = 2.5,
  disp.min = -2,
  disp.max = 2,
  assay = "ADT"
) +
  scale_fill_gradientn(
    limits = c(-2, 2),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Non-responders, proteins")
dev.off()

hd_vs_post_res_deg <- hd_vs_post_res_deg[hd_vs_post_res_deg$p_val_adj < 0.05, ]

hist(hd_post_resp_so[rownames(hd_vs_post_res_deg)[rownames(hd_vs_post_res_deg) %in% nahids_genes], ]@assays$RNA@scale.data)

hd_vs_post_res_deg <- hd_vs_post_res_deg[order(hd_vs_post_res_deg$avg_logFC), ]
rownames(hd_vs_post_res_deg)[rownames(hd_vs_post_res_deg) %in% nahids_genes]

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/5.rna_heatmap_hd_vs_post_res_deg.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  hd_post_resp_so,
  features = rownames(hd_vs_post_res_deg)[rownames(hd_vs_post_res_deg) %in% nahids_genes],
  # features = rownames(hd_vs_post_res_deg),
  group.by = 'type',
  size = 2.5,
  disp.min = -2,
  disp.max = 4
) +
  scale_fill_gradientn(
    limits = c(-3, 4),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Healthy Donors vs Post-Responders, RNA")
dev.off()

hd_vs_post_res_adt_deg <- hd_vs_post_res_adt_deg[hd_vs_post_res_adt_deg$p_val_adj < 0.05, ]

hist(hd_post_resp_so@assays$ADT@scale.data[rownames(hd_vs_post_res_adt_deg), ])

hd_vs_post_res_adt_deg <- hd_vs_post_res_adt_deg[order(hd_vs_post_res_adt_deg$avg_logFC), ]

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/6.adt_heatmap_hd_vs_post_res_deg.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  hd_post_resp_so,
  features = rownames(hd_vs_post_res_adt_deg),
  group.by = 'type',
  size = 2.5,
  disp.min = -3,
  disp.max = 3,
  assay = "ADT"
) +
  scale_fill_gradientn(
    limits = c(-3, 4),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Healthy Donors vs Post-Responders, proteins")
dev.off()

hd_vs_post_nres_deg <- hd_vs_post_nres_deg[hd_vs_post_nres_deg$p_val_adj < 0.05, ]

hist(hd_post_nresp_so[rownames(hd_vs_post_nres_deg)[rownames(hd_vs_post_nres_deg) %in% nahids_genes], ]@assays$RNA@scale.data)

hd_vs_post_nres_deg <- hd_vs_post_nres_deg[order(hd_vs_post_nres_deg$avg_logFC), ]
rownames(hd_vs_post_nres_deg)[rownames(hd_vs_post_nres_deg) %in% nahids_genes]

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/7.rna_heatmap_hd_vs_post_nres_deg.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_deg)[rownames(hd_vs_post_nres_deg) %in% nahids_genes],
  group.by = 'type',
  size = 2.5,
  disp.min = -2,
  disp.max = 3
) +
  scale_fill_gradientn(
    limits = c(-2, 3),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Healthy Donors vs Post-Non-responders, RNA")
dev.off()

hd_vs_post_nres_adt_deg <- hd_vs_post_nres_adt_deg[hd_vs_post_nres_adt_deg$p_val_adj < 0.05, ]

hist(hd_post_nresp_so@assays$ADT@scale.data[rownames(hd_vs_post_nres_adt_deg), ])

hd_vs_post_nres_adt_deg <- hd_vs_post_nres_adt_deg[order(hd_vs_post_nres_adt_deg$avg_logFC), ]

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/4.heatmaps/8.adt_heatmap_hd_vs_post_nres_deg.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DoHeatmap(
  hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_adt_deg),
  group.by = 'type',
  size = 2.5,
  disp.min = -3,
  disp.max = 3,
  assay = "ADT"
) +
  scale_fill_gradientn(
    limits = c(-3, 3),
    colors = colorRampPalette(c("midnightblue", "red3", "yellow"))(50),
    na.value = "white"
  ) +
  ggtitle("Healthy Donors vs Post-Non-responders, proteins")
dev.off()

# Violin graphs
vplots = VlnPlot(
  object = resp_so,
  features = rownames(res_pre_vs_post_deg)[rownames(res_pre_vs_post_deg) %in% nahids_genes],
  split.by = "time",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/1.violins_for_selected_genes_pre_vs_post_resp.pdf",
    sep = ""
  ),
  width = 16,
  height = 20
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Responders, RNA, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = resp_so,
  features = rownames(res_pre_vs_post_adt_deg),
  split.by = "time",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1,
  assay = "ADT"
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/2.violins_for_selected_genes_pre_vs_post_resp_adt.pdf",
    sep = ""
  ),
  width = 16,
  height = 11
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Responders, proteins, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = nresp_so,
  features = rownames(nres_pre_vs_post_deg)[rownames(nres_pre_vs_post_deg) %in% nahids_genes],
  split.by = "time",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/3.violins_for_selected_genes_pre_vs_post_nresp.pdf",
    sep = ""
  ),
  width = 16,
  height = 20
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Non-responders, RNA, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = nresp_so,
  features = rownames(nres_pre_vs_post_adt_deg),
  split.by = "time",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1,
  assay = "ADT"
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/4.violins_for_selected_genes_pre_vs_post_nresp_adt.pdf",
    sep = ""
  ),
  width = 16,
  height = 11
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Non-responders, proteins, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = hd_post_resp_so,
  features = rownames(hd_vs_post_res_deg)[rownames(hd_vs_post_res_deg) %in% nahids_genes],
  split.by = "type",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/5.violins_for_selected_genes_hd_vs_post_resp.pdf",
    sep = ""
  ),
  width = 16,
  height = 24
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Healthy donors vs post responders, RNA, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = hd_post_resp_so,
  features = rownames(hd_vs_post_res_adt_deg),
  split.by = "type",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1,
  assay = "ADT"
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/6.violins_for_selected_genes_hd_vs_post_resp_adt.pdf",
    sep = ""
  ),
  width = 16,
  height = 13
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Healthy donors vs post responders, proteins, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_deg)[rownames(hd_vs_post_nres_deg) %in% nahids_genes],
  split.by = "type",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/7.violins_for_selected_genes_hd_vs_post_nresp.pdf",
    sep = ""
  ),
  width = 16,
  height = 7
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Healthy donors vs post non-responders, RNA, SNN clusters', face = 'bold')
)
dev.off()

vplots = VlnPlot(
  object = hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_adt_deg),
  split.by = "type",
  split.plot = TRUE,
  combine = FALSE,
  pt.size = 0,
  ncol = 1,
  assay = "ADT"
)
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/5.violins/8.violins_for_selected_genes_hd_vs_post_nresp_adt.pdf",
    sep = ""
  ),
  width = 16,
  height = 20
)
cps <- CombinePlots(plots = vplots, ncol = 1)
annotate_figure(
  p = cps,
  top = ggpubr::text_grob(label = 'Healthy donors vs post non-responders, proteins, SNN clusters', face = 'bold')
)
dev.off()

# Feature plots
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/1.feature_plots_for_selected_genes_pre_vs_post_resp.pdf",
    sep = ""
  ),
  width = 12,
  height = 20
)
FeaturePlot(
  resp_so,
  features = rownames(res_pre_vs_post_deg)[rownames(res_pre_vs_post_deg) %in% nahids_genes],
  split.by = "time",
  max.cutoff = 3
)
dev.off()

DefaultAssay(resp_so) <- "ADT"
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/2.feature_plots_for_selected_genes_pre_vs_post_resp_adt.pdf",
    sep = ""
  ),
  width = 12,
  height = 12
)
FeaturePlot(
  resp_so,
  features = rownames(res_pre_vs_post_adt_deg),
  split.by = "time",
  max.cutoff = 3
)
dev.off()
DefaultAssay(resp_so) <- "RNA"

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/3.feature_plots_for_selected_genes_pre_vs_post_nresp.pdf",
    sep = ""
  ),
  width = 12,
  height = 36
)
FeaturePlot(
  nresp_so,
  features = rownames(nres_pre_vs_post_deg)[rownames(nres_pre_vs_post_deg) %in% nahids_genes],
  split.by = "time",
  max.cutoff = 3
)
dev.off()

DefaultAssay(nresp_so) <- "ADT"
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/4.feature_plots_for_selected_genes_pre_vs_post_nresp_adt.pdf",
    sep = ""
  ),
  width = 12,
  height = 12
)
FeaturePlot(
  nresp_so,
  features = rownames(nres_pre_vs_post_adt_deg),
  split.by = "time",
  max.cutoff = 3
)
dev.off()
DefaultAssay(nresp_so) <- "RNA"

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/5.feature_plots_for_selected_genes_hd_vs_post_res.pdf",
    sep = ""
  ),
  width = 12,
  height = 28
)
FeaturePlot(
  hd_post_resp_so,
  features = rownames(hd_vs_post_res_deg)[rownames(hd_vs_post_res_deg) %in% nahids_genes],
  split.by = "type",
  max.cutoff = 3
)
dev.off()

DefaultAssay(hd_post_resp_so) <- "ADT"
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/6.feature_plots_for_selected_genes_hd_vs_post_res_adt.pdf",
    sep = ""
  ),
  width = 12,
  height = 16
)
FeaturePlot(
  hd_post_resp_so,
  features = rownames(hd_vs_post_res_adt_deg),
  split.by = "type",
  max.cutoff = 3
)
dev.off()
DefaultAssay(hd_post_resp_so) <- "RNA"

pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/7.feature_plots_for_selected_genes_hd_vs_post_nres.pdf",
    sep = ""
  ),
  width = 12,
  height = 8
)
FeaturePlot(
  hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_deg)[rownames(hd_vs_post_nres_deg) %in% nahids_genes],
  split.by = "type",
  max.cutoff = 3
)
dev.off()

DefaultAssay(hd_post_nresp_so) <- "ADT"
pdf(
  paste(
    docs_dir,
    "plots/b_corrected/6.feature_plots/8.feature_plots_for_selected_genes_hd_vs_post_nres_adt.pdf",
    sep = ""
  ),
  width = 12,
  height = 24
)
FeaturePlot(
  hd_post_nresp_so,
  features = rownames(hd_vs_post_nres_adt_deg),
  split.by = "type",
  max.cutoff = 3
)
dev.off()
DefaultAssay(hd_post_nresp_so) <- "RNA"
