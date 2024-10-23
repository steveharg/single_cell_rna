# - performs GSVA enrichment when comparing two cell populations
require(Seurat)
require(GSVA)
require(GSEABase)
require(GSVAdata)
data(c2BroadSets) #http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2
require(Biobase)
require(genefilter)
require(limma)
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(corrplot)
require(biomaRt)
require(config)
source('gsva_heatmap.r')

# Initial config
data_dir = "../data/"
docs_dir = "../docs/"

load('../data/umap_pre_enrichment.rdata')

#################################
# Perform GSVA analysis
# compare pathways identified in each time-point against each other
# pathways identified from statistically significantly differentially expressed genes
#################################

#################################################
# Pre vs Post

# KEGG pathways
# c2BroadSets <- c2BroadSets[which(names(c2BroadSets) %in% names(c2BroadSets)[grep("KEGG_|BIOCARTA_|REACTOME_", names(c2BroadSets))])]
counts_pre_post_so = subset(counts_combined_so, subset = (time == 'post' |
                                                            time == 'pre'))
topMatrix <- counts_pre_post_so@assays$integrated@scale.data

#Convert the HGNC names to Entrez IDs (eliminate non-matches where necessary)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annots <- getBM(
  mart = mart,
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filter = "hgnc_symbol",
  values = rownames(topMatrix),
  uniqueRows = TRUE
)
annots <- annots[!duplicated(annots[, 1]), ]
topMatrix <- topMatrix[which(rownames(topMatrix) %in% annots[, 1]), ]
annots <- annots[which(annots[, 1] %in% rownames(topMatrix)), ]
topMatrix <- topMatrix[match(annots[, 1], rownames(topMatrix)), ]
rownames(topMatrix) <- annots[, 2]

#Perform GSVA
#topMatrixGSVA <- gsva(topMatrix, c2BroadSets, method = "zscore", min.sz=5, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)
topMatrixGSVA <- gsva(
  topMatrix,
  c2BroadSets,
  method = "gsva",
  min.sz = 10,
  max.sz = 999999,
  abs.ranking = FALSE,
  verbose = TRUE
)

metadata = data.frame(counts_pre_post_so$time)
colnames(metadata) = "Time"

topMatrixGSVAsubset = t(topMatrixGSVA[, colnames(counts_pre_post_so)])

# Verify the order of the cell IDs match in both metadata and topMatrixGSVA
all_match = TRUE
for (i in 1:dim(metadata)[1]) {
  if (rownames(metadata)[i] != rownames(topMatrixGSVAsubset)[i]) {
    all_match = FALSE
  }
}
if (!all_match) {
  stop('The order of the cell IDs in both metadata and topMatrixGSVA does not match')
}

modeling <- data.frame(metadata, topMatrixGSVAsubset)
design <- model.matrix( ~ modeling$Time)
colnames(design) <- c("Intercept", "Pre vs Post")
fit <- lmFit(t(modeling[, 2:ncol(modeling)]), design)
fit <- eBayes(fit)
pre_vs_post <- topTable(
  fit,
  coef = "Pre vs Post",
  number = Inf,
  p.value = 1,
  adjust = "BH"
)

out = data.frame(pathway = rownames(pre_vs_post), pre_vs_post[, 1:ncol(pre_vs_post)])
write.table(
  x = out,
  file = paste0(data_dir, "pre_vs_post_ranked_pathways.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

colours <- list("Time" = c("Pre" = "royalblue", "Post" = "red2"))

hm_results = gsva_heatmap(
  topRankedPathways = pre_vs_post,
  seurat_sample = counts_pre_post_so,
  metadata_column = metadata$Time,
  colours
)

pdf(
  paste0(
    docs_dir,
    "plots/integrated/pathways/Pre_vs_Post.ComplexHeatmap.pdf"
  ),
  width = 16,
  height = 8
)
draw(
  hm_results[[1]] + hm_results[[2]],
  heatmap_legend_side = "bottom",
  annotation_legend_side = "top",
  row_sub_title_side = "left"
)
dev.off()

p <- aggregate(modeling[, 2:ncol(modeling)], by = modeling[1], mean)
rownames(p) <- p[, 1]
p <- p[, -1]
filt <- rownames(pre_vs_post[which(pre_vs_post$adj.P.Val <= 0.01), ])
p <- p[, filt]
p <- t(p) / max(abs(range(p))) # scale between -1 and +1
pick.col <- brewer.pal(10, "RdBu")
col <- c(colorRampPalette(rev(pick.col))(200))
pdf(
  paste0(docs_dir, "plots/integrated/pathways/pre_vs_post_GSVA.pdf"),
  width = 30,
  height = 5
)
par(xpd = TRUE)
corrplot(
  t(p),
  method = "circle",
  order = "original",
  addgrid.col = "grey60",
  tl.col = "black",
  col = col,
  cl.cex = 1.0,
  cl.pos = "b",
  cl.ratio = 1.5,
  cl.lim = c(-1, 1),
  cl.length = 5,
  tl.cex = 0.6,
  tl.srt = 35,
  mar = c(0.1, 0.1, 0.1, 20)
)
legend(
  "top",
  bty = "n",
  cex = 1.0,
  title = "Pathway activation (adj p<0.01)",
  c("High", "Mean", "Low"),
  fill = c(col[length(col) - 10], col[length(col) / 2], col[11]),
  horiz = TRUE
)
dev.off()

#################################################
# HD vs Pre

# KEGG pathways
# c2BroadSets <- c2BroadSets[which(names(c2BroadSets) %in% names(c2BroadSets)[grep("KEGG_|BIOCARTA_|REACTOME_", names(c2BroadSets))])]
counts_hd_pre_so = subset(counts_combined_so,
                          subset = (hd_or_pre == 'hd' | hd_or_pre == 'pre'))
topMatrix <- counts_hd_pre_so@assays$integrated@scale.data

#Convert the HGNC names to Entrez IDs (eliminate non-matches where necessary)
annots <- getBM(
  mart = mart,
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filter = "hgnc_symbol",
  values = rownames(topMatrix),
  uniqueRows = TRUE
)
annots <- annots[!duplicated(annots[, 1]), ]
topMatrix <- topMatrix[which(rownames(topMatrix) %in% annots[, 1]), ]
annots <- annots[which(annots[, 1] %in% rownames(topMatrix)), ]
topMatrix <- topMatrix[match(annots[, 1], rownames(topMatrix)), ]
rownames(topMatrix) <- annots[, 2]

#Perform GSVA
#topMatrixGSVA <- gsva(topMatrix, c2BroadSets, method = "zscore", min.sz=5, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)
topMatrixGSVA <- gsva(
  topMatrix,
  c2BroadSets,
  method = "gsva",
  min.sz = 10,
  max.sz = 999999,
  abs.ranking = FALSE,
  verbose = TRUE
)

metadata = data.frame(counts_hd_pre_so$hd_or_pre)
colnames(metadata) = "hd_or_pre"

topMatrixGSVAsubset = t(topMatrixGSVA[, colnames(counts_hd_pre_so)])

# Verify the order of the cell IDs match in both metadata and topMatrixGSVA
all_match = TRUE
for (i in 1:dim(metadata)[1]) {
  if (rownames(metadata)[i] != rownames(topMatrixGSVAsubset)[i]) {
    all_match = FALSE
  }
}
if (!all_match) {
  stop('The order of the cell IDs in both metadata and topMatrixGSVA does not match')
}

modeling <- data.frame(metadata, topMatrixGSVAsubset)
design <- model.matrix( ~ modeling$hd_or_pre)
colnames(design) <- c("Intercept", "HD vs Pre")
fit <- lmFit(t(modeling[, 2:ncol(modeling)]), design)
fit <- eBayes(fit)
hd_vs_pre <- topTable(
  fit,
  coef = "HD vs Pre",
  number = Inf,
  p.value = 1,
  adjust = "BH"
)

out = data.frame(pathway = rownames(hd_vs_pre), hd_vs_pre[, 1:ncol(hd_vs_pre)])
write.table(
  x = out,
  file = paste0(data_dir, "hd_vs_pre_ranked_pathways.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

colours <- list("Time" = c("HD" = "royalblue", "Pre" = "red2"))

hm_results = gsva_heatmap(
  topRankedPathways = hd_vs_pre,
  seurat_sample = counts_hd_pre_so,
  metadata_column = metadata$hd_or_pre,
  colours
)

pdf(
  paste0(
    docs_dir,
    "plots/integrated/pathways/hd_vs_pre.ComplexHeatmap.pdf"
  ),
  width = 16,
  height = 8
)
draw(
  hm_results[[1]] + hm_results[[2]],
  heatmap_legend_side = "bottom",
  annotation_legend_side = "top",
  row_sub_title_side = "left"
)
dev.off()

p <- aggregate(modeling[, 2:ncol(modeling)], by = modeling[1], mean)
rownames(p) <- p[, 1]
p <- p[, -1]
filt <- rownames(hd_vs_pre[which(hd_vs_pre$adj.P.Val <= 0.01), ])
p <- p[, filt]
p <- t(p) / max(abs(range(p))) # scale between -1 and +1
pick.col <- brewer.pal(10, "RdBu")
col <- c(colorRampPalette(rev(pick.col))(200))
pdf(
  paste0(docs_dir, "plots/integrated/pathways/hd_vs_pre_GSVA.pdf"),
  width = 30,
  height = 5
)
par(xpd = TRUE)
corrplot(
  t(p),
  method = "circle",
  order = "original",
  addgrid.col = "grey60",
  tl.col = "black",
  col = col,
  cl.cex = 1.0,
  cl.pos = "b",
  cl.ratio = 1.5,
  cl.lim = c(-1, 1),
  cl.length = 5,
  tl.cex = 0.6,
  tl.srt = 35,
  mar = c(0.1, 0.1, 0.1, 20)
)
legend(
  "top",
  bty = "n",
  cex = 1.0,
  title = "Pathway activation (adj p<0.01)",
  c("High", "Mean", "Low"),
  fill = c(col[length(col) - 10], col[length(col) / 2], col[11]),
  horiz = TRUE
)
dev.off()

#################################################
# HD vs Post

# KEGG pathways
# c2BroadSets <- c2BroadSets[which(names(c2BroadSets) %in% names(c2BroadSets)[grep("KEGG_|BIOCARTA_|REACTOME_", names(c2BroadSets))])]
counts_hd_post_so = subset(counts_combined_so,
                           subset = (hd_or_post == 'hd' | hd_or_post == 'post'))
topMatrix <- counts_hd_post_so@assays$integrated@scale.data

#Convert the HGNC names to Entrez IDs (eliminate non-matches where necessary)
annots <- getBM(
  mart = mart,
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filter = "hgnc_symbol",
  values = rownames(topMatrix),
  uniqueRows = TRUE
)
annots <- annots[!duplicated(annots[, 1]), ]
topMatrix <- topMatrix[which(rownames(topMatrix) %in% annots[, 1]), ]
annots <- annots[which(annots[, 1] %in% rownames(topMatrix)), ]
topMatrix <- topMatrix[match(annots[, 1], rownames(topMatrix)), ]
rownames(topMatrix) <- annots[, 2]

#Perform GSVA
#topMatrixGSVA <- gsva(topMatrix, c2BroadSets, method = "zscore", min.sz=5, max.sz=999999, abs.ranking=FALSE, verbose=TRUE)
topMatrixGSVA <- gsva(
  topMatrix,
  c2BroadSets,
  method = "gsva",
  min.sz = 10,
  max.sz = 999999,
  abs.ranking = FALSE,
  verbose = TRUE
)

metadata = data.frame(counts_hd_post_so$hd_or_post)
colnames(metadata) = "hd_or_post"

topMatrixGSVAsubset = t(topMatrixGSVA[, colnames(counts_hd_post_so)])

# Verify the order of the cell IDs match in both metadata and topMatrixGSVA
all_match = TRUE
for (i in 1:dim(metadata)[1]) {
  if (rownames(metadata)[i] != rownames(topMatrixGSVAsubset)[i]) {
    all_match = FALSE
  }
}
if (!all_match) {
  stop('The order of the cell IDs in both metadata and topMatrixGSVA does not match')
}

modeling <- data.frame(metadata, topMatrixGSVAsubset)
design <- model.matrix( ~ modeling$hd_or_post)
colnames(design) <- c("Intercept", "HD vs Post")
fit <- lmFit(t(modeling[, 2:ncol(modeling)]), design)
fit <- eBayes(fit)
hd_vs_post <- topTable(
  fit,
  coef = "HD vs Post",
  number = Inf,
  p.value = 1,
  adjust = "BH"
)

out = data.frame(pathway = rownames(hd_vs_post), hd_vs_post[, 1:ncol(hd_vs_post)])
write.table(
  x = out,
  file = paste0(data_dir, "hd_vs_post_ranked_pathways.tsv"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

colours <- list("Time" = c("HD" = "royalblue", "Post" = "red2"))

hm_results = gsva_heatmap(
  topRankedPathways = hd_vs_post,
  seurat_sample = counts_hd_post_so,
  metadata_column = metadata$hd_vs_post,
  colours
)

pdf(
  paste0(
    docs_dir,
    "plots/integrated/pathways/hd_vs_post.ComplexHeatmap.pdf"
  ),
  width = 16,
  height = 8
)
draw(
  hm_results[[1]] + hm_results[[2]],
  heatmap_legend_side = "bottom",
  annotation_legend_side = "top",
  row_sub_title_side = "left"
)
dev.off()

p <- aggregate(modeling[, 2:ncol(modeling)], by = modeling[1], mean)
rownames(p) <- p[, 1]
p <- p[, -1]
filt <- rownames(hd_vs_post[which(hd_vs_post$adj.P.Val <= 0.01), ])
p <- p[, filt]
p <- t(p) / max(abs(range(p))) # scale between -1 and +1
pick.col <- brewer.pal(10, "RdBu")
col <- c(colorRampPalette(rev(pick.col))(200))
pdf(
  paste0(docs_dir, "plots/integrated/pathways/hd_vs_post_GSVA.pdf"),
  width = 30,
  height = 5
)
par(xpd = TRUE)
corrplot(
  t(p),
  method = "circle",
  order = "original",
  addgrid.col = "grey60",
  tl.col = "black",
  col = col,
  cl.cex = 1.0,
  cl.pos = "b",
  cl.ratio = 1.5,
  cl.lim = c(-1, 1),
  cl.length = 5,
  tl.cex = 0.6,
  tl.srt = 35,
  mar = c(0.1, 0.1, 0.1, 20)
)
legend(
  "top",
  bty = "n",
  cex = 1.0,
  title = "Pathway activation (adj p<0.01)",
  c("High", "Mean", "Low"),
  fill = c(col[length(col) - 10], col[length(col) / 2], col[11]),
  horiz = TRUE
)
dev.off()
