library(garnett)
library(org.Hs.eg.db)
library(Seurat)
library(cowplot)

rm(list = ls())

load('../data/processed_cts.rdata')

pbmc_classifier = readRDS("../data/hsPBMC.RDS")

# mat <- Matrix::readMM(system.file("extdata", "exprs_sparse.mtx", package = "garnett"))
# fdata <- read.table(system.file("extdata", "fdata.txt", package = "garnett"))
# pdata <- read.table(system.file("extdata", "pdata.txt", package = "garnett"),
#                     sep="\t")
# row.names(mat) <- row.names(fdata)
# colnames(mat) <- row.names(pdata)
#
# # create a new CDS object
# pd <- new("AnnotatedDataFrame", data = pdata)
# fd <- new("AnnotatedDataFrame", data = fdata)
# pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
#                            phenoData = pd,
#                            featureData = fd)
#
# # generate size factors for normalization later
# pbmc_cds <- estimateSizeFactors(pbmc_cds)

fdata2 = data.frame(gene_short_name = rownames(counts_combined_so@assays$integrated@scale.data))
rownames(fdata2) = fdata2$gene_short_name
fd2 <- new("AnnotatedDataFrame", data = fdata2)

pdata2 = counts_combined_so[["rnaClusterID"]]
pdata2$garnett_cluster = pdata2$rnaClusterID
pdata2 = pdata2[-1]
pd2 <- new("AnnotatedDataFrame", data = pdata2)

counts_combined_cds = newCellDataSet(
  counts_combined_so@assays$integrated@scale.data,
  phenoData = pd2,
  featureData = fd2
)


counts_combined_cds = estimateSizeFactors(counts_combined_cds)

counts_combined_cds <- classify_cells(
  counts_combined_cds,
  pbmc_classifier,
  db = org.Hs.eg.db,
  cluster_extend = TRUE,
  cds_gene_id_type = "SYMBOL"
)


head(pData(counts_combined_cds))

counts_combined_so <- AddMetaData(object = counts_combined_so,
                                  metadata = counts_combined_cds$cell_type,
                                  col.name = 'cell_type')

counts_combined_so <- AddMetaData(
  object = counts_combined_so,
  metadata = counts_combined_cds$cluster_ext_type,
  col.name = 'cluster_ext_type'
)

p1 = DimPlot(counts_combined_so,
             reduction = "umap",
             group.by = "patient") + ggtitle("All projects merged, patients")
p2 = DimPlot(counts_combined_so,
             reduction = "umap",
             group.by = "time") + ggtitle("All projects merged, pre vs post")
p3 = DimPlot(counts_combined_so,
             reduction = "umap",
             group.by = "stimulation") + ggtitle("All projects merged, stimulation")
p4 = DimPlot(
  counts_combined_so,
  reduction = "umap",
  group.by = "rnaClusterID",
  label = TRUE,
  do.return = TRUE
) + ggtitle("All projects merged, clusters")
p5 = DimPlot(
  counts_combined_so,
  reduction = "umap",
  group.by = "cluster_ext_type",
  label = TRUE,
  do.return = TRUE
) + ggtitle("All projects merged, extended cell types")
pdf(
  paste(
    docs_dir,
    "plots/integrated/analyses/umap_with_cell_types.pdf",
    sep = ""
  ),
  width = 16,
  height = 8
)
plot_grid(p1, p2, p3, p4, p5)
dev.off()

uniq_clusters_cell_types = unique(
  data.frame(
    cluster = counts_combined_cds$garnett_cluster,
    cell_type = counts_combined_cds$cell_type
  )
)
uniq_clusters_ext_types = unique(
  data.frame(
    cluster = counts_combined_cds$garnett_cluster,
    ext_type = counts_combined_cds$cluster_ext_type
  )
)

uniq_clusters_cell_types[order(uniq_clusters_cell_types$cluster), ]
uniq_clusters_ext_types[order(uniq_clusters_ext_types$cluster), ]