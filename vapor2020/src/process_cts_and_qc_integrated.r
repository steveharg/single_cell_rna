# - a seurat process tailored for Nahid's 12 VAPOR RNA + protein data

# Clears workspace
rm(list = ls())

require("biomaRt")
require(doParallel)
library(reshape2)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
require("gplots")
source("shorten_gene_abseq_names.r")

# Initial config
rna_clusters_resolution <- 0.8
# adt_clusters_resolution <- 0.2
adt_clusters_resolution <- 0.8
counts_file_header <- TRUE
counts_file_skip <- 7
sample_tag_calls_header <- TRUE
sample_tag_calls_skip <- 7
data_dir <- "../data/"
data_in_dir <- "../data/input/"
docs_dir <- "../docs/"

CreateProject <- setRefClass(
  "Project",
  fields = list(
    name = "character",
    counts_file = "character",
    counts = "matrix",
    samples = "list",
    so = "Seurat"
  )
)

CreateSample <- setRefClass(
  "Sample",
  fields = list(
    sample_tag = "character",
    sample_tag_calls_file = "character",
    cells = "character",
    so = "Seurat"
  )
)

samples_meta <- read.table(
  file = paste0(data_in_dir, "sample_metadata.tsv"),
  sep = "\t",
  header = TRUE
)

projects <- list()
proj_idx <- 1
short_patient_gene_abseq_names <- NULL

for (p in  unique(samples_meta$project)) {
  proj_samples <- which(samples_meta$project == p)
  
  project <- CreateProject(name = p, counts_file = as.character(unique(samples_meta[which(samples_meta$project == p), "counts_file"])))
  
  counts <- read.csv(paste0(data_in_dir, project$counts_file),
                     header = counts_file_header,
                     skip = counts_file_skip)
  rownames(counts) <- counts[, 1]
  counts <- counts[, -1]
  counts <- t(counts)
  if (is.null(short_patient_gene_abseq_names)) {
    short_patient_gene_abseq_names <- as.character(shorten_gene_abseq_names(counts))
  }
  rownames(counts) <- short_patient_gene_abseq_names
  
  project$counts <- counts
  # Create seurat object (separate assays for rna and proteins)
  project$so <- CreateSeuratObject(counts = project$counts[grep('pAbO', invert = TRUE, rownames(project$counts)), ])
  project$so[["ADT"]] <- CreateAssayObject(counts = project$counts[grep('pAbO', rownames(project$counts)), ])
  
  i <- 1
  samples <- list()
  sample_tag_calls = list()
  
  for (proj_sample in proj_samples) {
    sample <- CreateSample(
      sample_tag = as.character(samples_meta[proj_sample, "sample_tag"]),
      sample_tag_calls_file = as.character(samples_meta[proj_sample, "sample_tag_calls_file"])
    )
    
    st_file <- sample$sample_tag_calls_file
    
    if (!st_file %in% names(sample_tag_calls)) {
      sample_tag_calls[[st_file]] <- read.csv(
        header = sample_tag_calls_header,
        file = paste0(data_in_dir, st_file),
        skip = sample_tag_calls_skip
      )
      rownames(sample_tag_calls[[st_file]]) <- sample_tag_calls[[st_file]][, 1]
      sample_tag_calls[[st_file]] <- sample_tag_calls[[st_file]][, -1]
    }
    
    sample$cells <- colnames(project$counts)[which(sample_tag_calls[[st_file]]$Sample_Tag == sample$sample_tag)]
    
    if (length(sample$cells) > 0) {
      # Metadata allowed types
      # $stimulation = "stim" / "unstim"
      # $time = "pre" / "post" / "na"
      # $type = "resp" / "nonresp" / "hd"
      
      # Get subsets for each sample, both RNA and ADT
      sample$so <- subset(x = project$so, cells = sample$cells)
      sample$so$stimulation <- as.character(samples_meta[proj_sample, "stimulation"])
      sample$so$time <- as.character(samples_meta[proj_sample, "time"])
      sample$so$patient <- as.character(samples_meta[proj_sample, "patient"])
      sample$so$type <- as.character(samples_meta[proj_sample, "type"])
      sample$so$project <- as.character(samples_meta[proj_sample, "project"])
      sample$so <- NormalizeData(sample$so)
      sample$so <- NormalizeData(sample$so,
                                 assay = "ADT",
                                 normalization.method = "CLR")
      sample$so <- FindVariableFeatures(sample$so)
      
      samples[[i]] <- sample
    }
    
    i <- i + 1
  }
  
  project$samples <- samples
  projects[[proj_idx]] <- project
  proj_idx <- proj_idx + 1
}

# Check gene lists match between projects, and save
# full gene/AbSeq names before shortening
all_genes_match <- TRUE
ensmbl_patient_gene_abseq_names <- rownames(projects[[1]]$counts)
for (project in projects) {
  if (!all(rownames(project$counts) == ensmbl_patient_gene_abseq_names)) {
    all_genes_match <- FALSE
  }
}
if (!all_genes_match) {
  stop("gene names in counts file do not match")
}

all_sample_sos <- list()
for (project in projects) {
  for (sample in project$samples) {
    all_sample_sos <- c(all_sample_sos, sample$so)
  }
}

# More advanced merge - does this work for Nahid's dataset ?
# Integrate samples
rna_anchors <- FindIntegrationAnchors(object.list = all_sample_sos[c(1, 3:length(all_sample_sos))], k.filter = 100)
rna_integrated_so <- IntegrateData(anchorset = rna_anchors, dims = 1:30)

all_sample_adt_sos = list()
for (sample_so in all_sample_sos) {
  DefaultAssay(sample_so) <- "ADT"
  all_sample_adt_sos <- c(all_sample_adt_sos, sample_so)
}

adt_anchors <- FindIntegrationAnchors(object.list = all_sample_adt_sos[c(1, 3:length(all_sample_adt_sos))], k.filter = 100)
adt_integrated_so <- IntegrateData(anchorset = adt_anchors, dims = 1:8)

##### ******** Should we use merge.data = FALSE ? that would mean don't normalise before, do it after the merge
merged_so = merge(all_sample_sos[[1]], all_sample_sos[2:length(all_sample_sos)])

save.image('../data/processed_cts.rdata')

##### UP TO HERE ####

# Examine histograms of RNA and protein - RNA should be approx negative binomial, protein should be something different
# (approx normal or binomial)
pdf(
  paste(docs_dir, "plots/integrated/qc/rna_hist.pdf", sep = ""),
  width = 16,
  height = 8
)
hist(as.matrix(rna_integrated_so@assays$integrated@data),
     breaks = 100)
dev.off()
pdf(
  paste(docs_dir, "plots/integrated/qc/adt_hist.pdf", sep = ""),
  width = 16,
  height = 8
)
hist(as.matrix(adt_integrated_so@assays$integrated@data),
     breaks = 20)
dev.off()

# Quality checks

# Visualize QC metrics as a violin plot
# Check that number of distinct genes per cell ('nFeature_RNA') and
# the sum of the counts for a cell ('nCount_RNA') look reasonable.
pdf(
  paste(docs_dir, "plots/integrated/qc/rna_qc_violins.pdf", sep = ""),
  width = 16,
  height = 8
)
VlnPlot(merged_so,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 3)
dev.off()

pdf(
  paste(docs_dir, "plots/integrated/qc/adt_qc_violins.pdf", sep = ""),
  width = 16,
  height = 8
)
VlnPlot(
  merged_so,
  features = c("nFeature_ADT", "nCount_ADT"),
  ncol = 3,
  assay = "ADT"
)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

pdf(
  paste(
    docs_dir,
    "plots/integrated/qc/rna_qc_feature_scatters.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
FeatureScatter(merged_so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

pdf(
  paste(
    docs_dir,
    "plots/integrated/qc/adt_qc_feature_scatters.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
FeatureScatter(merged_so, feature1 = "nCount_ADT", feature2 = "nFeature_ADT")
dev.off()

# choose ~1k variable features
merged_so <- FindVariableFeatures(merged_so)
rna_integrated_so <- FindVariableFeatures(rna_integrated_so)

# Identify the 10 most highly variable genes
top10_rna_merged <- head(VariableFeatures(merged_so), 10)
top10_rna_integrated <- head(VariableFeatures(rna_integrated_so), 10)

# standard scaling (no regression)
merged_so <- ScaleData(merged_so)
rna_integrated_so <- ScaleData(rna_integrated_so)

# Run PCA, select PCs for tSNE / UMAP visualization and graph-based clustering
merged_so <- RunPCA(merged_so, verbose = FALSE)
rna_integrated_so <- RunPCA(rna_integrated_so, verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(merged_so[["pca"]], dims = 1:5, nfeatures = 5)
print(rna_integrated_so[["pca"]], dims = 1:5, nfeatures = 5)

pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca_loadings.pdf", sep = ""),
  width = 16,
  height = 8
)
VizDimLoadings(merged_so, dims = 1:2, reduction = "pca")
dev.off()

pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca.pdf", sep = ""),
  width = 16,
  height = 8
)
DimPlot(merged_so, reduction = "pca")
dev.off()
pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca1_heatmap.pdf", sep = ""),
  width = 16,
  height = 8
)
DimHeatmap(merged_so,
           dims = 1,
           cells = 500,
           balanced = TRUE)
dev.off()
pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca_heatmap.pdf", sep = ""),
  width = 16,
  height = 8
)
DimHeatmap(merged_so,
           dims = 1:15,
           cells = 500,
           balanced = TRUE)
dev.off()

pdf(
  paste(
    docs_dir,
    "plots/integrated/qc/rna_pca_by_stimulation.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
DimPlot(merged_so, reduction = "pca", group.by = "stimulation")
dev.off()

pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca_by_time.pdf", sep = ""),
  width = 16,
  height = 8
)
DimPlot(merged_so, reduction = "pca", group.by = "time")
dev.off()

pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca_by_patient.pdf", sep = ""),
  width = 16,
  height = 8
)
DimPlot(merged_so, reduction = "pca", group.by = "patient")
dev.off()

pdf(
  paste(docs_dir, "plots/integrated/qc/rna_pca_elbow.pdf", sep = ""),
  width = 16,
  height = 8
)
ElbowPlot(merged_so, ndims = 50)
dev.off()
# based on the elbow plot, select num_pcs PCs for UMAP visualization and graph-based clustering
num_pcs = 20

merged_so <- FindNeighbors(merged_so, dims = 1:num_pcs)
merged_so <- FindClusters(merged_so, resolution = rna_clusters_resolution)
merged_so <- RunUMAP(merged_so, reduction = "pca", dims = 1:num_pcs)
# merged_so <- RunTSNE(merged_so, reduction = "pca", dims = 1:num_pcs)
rna_integrated_so <- FindNeighbors(rna_integrated_so, dims = 1:num_pcs)
rna_integrated_so <- FindClusters(rna_integrated_so, resolution = rna_clusters_resolution)
rna_integrated_so <- RunUMAP(rna_integrated_so,
                             reduction = "pca",
                             dims = 1:num_pcs)
# rna_integrated_so <- RunTSNE(rna_integrated_so, reduction = "pca", dims = 1:num_pcs)

# stim_cells = Cells(merged_so[,which(merged_so@meta.data$stimulation == "stim")])
# unstim_cells = Cells(merged_so[,which(merged_so@meta.data$stimulation == "unstim")])
# pre_cells = Cells(merged_so[,which(merged_so@meta.data$time == "pre")])
# post_cells = Cells(merged_so[,which(merged_so@meta.data$time == "post")])
#
# patient18_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "18")])
# patient19_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "19")])
# patient24_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "24")])
# hd_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "hd")])
# patient20_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "20")])
# patient21_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "21")])
# patient22_cells = Cells(merged_so[,which(merged_so@meta.data$patient == "22")])

# # Now we can repeat the preprocessing (normalization and scaling) steps that we typically run
# # with RNA, but modifying the 'assay' argument.  For CITE-seq data, we do not recommend typical
# # LogNormalization. Instead, we use a centered log-ratio (CLR) normalization, computed
# # independently for each feature.  This is a slightly improved procedure from the original
# # publication, and we will release more advanced versions of CITE-seq normalizations soon.
# merged_so <- NormalizeData(merged_so, assay = "ADT", normalization.method = "CLR")


merged_so <- ScaleData(merged_so, assay = "ADT")
adt_integrated_so <- ScaleData(adt_integrated_so)

hist(merged_so@assays$RNA@scale.data)
hist(merged_so@assays$ADT@scale.data)
hist(rna_integrated_so@assays$integrated@scale.data)
hist(adt_integrated_so@assays$integrated@scale.data)


num_clusters_merged = length(levels(merged_so@active.ident))
cluster_cells_merged <- list()
for (i in 0:(num_clusters_merged - 1)) {
  cluster_cells_merged[[i + 1]] <- which(merged_so@active.ident == i)
}

markers_per_cluster_merged = list()
for (cluster in 0:(num_clusters_merged - 1)) {
  markers_per_cluster_merged[[cluster + 1]] = FindMarkers(object = merged_so, ident.1 = cluster)
}

all_markers_merged <- FindAllMarkers(merged_so) # are these the same as markers_per_cluster_merged ?

num_clusters_integrated = length(levels(rna_integrated_so@active.ident))
cluster_cells_integrated <- list()
for (i in 0:(num_clusters_integrated - 1)) {
  cluster_cells_integrated[[i + 1]] <- which(rna_integrated_so@active.ident == i)
}

markers_per_cluster_integrated = list()
for (cluster in 0:(num_clusters_integrated - 1)) {
  markers_per_cluster_integrated[[cluster + 1]] = FindMarkers(object = rna_integrated_so, ident.1 = cluster)
}

all_markers_integrated <- FindAllMarkers(rna_integrated_so) # are these the same as markers_per_cluster_integrated ?

# Because we're going to be working with the ADT data extensively, we're going to switch the
# default assay to the 'CITE' assay.  This will cause all functions to use ADT data by default,
# rather than requiring us to specify it each time
DefaultAssay(merged_so) <- "ADT"

# # based on the elbow plot, select num_pcs PCs for UMAP visualization and graph-based clustering
# num_pcs = 6

# Since we only have 10 markers, instead of doing PCA, we'll just use a standard euclidean
# distance matrix here.  Also, this provides a good opportunity to demonstrate how to do
# visualization and clustering using a custom distance matrix in Seurat
adt.data_merged <- GetAssayData(merged_so, slot = "data")
adt.dist_merged <- dist(t(adt.data_merged))
adt.data_integrated <- as.matrix(GetAssayData(adt_integrated_so, slot = "data"))
adt.dist_integrated <- dist(t(adt.data_integrated))

# Before we recluster the data on ADT levels, we'll stash the RNA cluster IDs for later
merged_so[["rnaClusterID"]] <- Idents(merged_so)
adt_integrated_so[["rnaClusterID"]] <- Idents(rna_integrated_so)

# Now, we rerun UMAP using our distance matrix defined only on ADT (protein) levels.
merged_so[["umap_adt"]] <- RunUMAP(adt.dist_merged, assay = "ADT", reduction.key = "adtUMAP_")
adt_integrated_so[["umap_adt"]] <- RunUMAP(adt.dist_integrated,
                                           assay = "integrated",
                                           reduction.key = "adtUMAP_")

# Fix strange bug - UMAP doesn't assign rownames, even though TSNE does
rownames(merged_so[['umap_adt']]@cell.embeddings) <- colnames(merged_so)
rownames(adt_integrated_so[['umap_adt']]@cell.embeddings) <- colnames(adt_integrated_so)

merged_so[["adt_snn"]] <- FindNeighbors(adt.dist_merged)$snn
adt_integrated_so[["adt_snn"]] <- FindNeighbors(adt.dist_integrated)$snn
merged_so <- FindClusters(merged_so, resolution = adt_clusters_resolution, graph.name = "adt_snn")
adt_integrated_so <- FindClusters(adt_integrated_so,
                                  resolution = adt_clusters_resolution,
                                  graph.name = "adt_snn")

# Reset default assays
DefaultAssay(merged_so) <- "RNA"

merged_so@meta.data$time = factor(merged_so@meta.data$time, levels = c("pre", "post"))
rna_integrated_so@meta.data$time = factor(rna_integrated_so@meta.data$time, levels = c("pre", "post"))
adt_integrated_so@meta.data$time = factor(rna_integrated_so@meta.data$time, levels = c("pre", "post"))

save.image('../data/processed_cts.rdata')
