#################################
# Enrichment via topGO and KEGGprofile
#################################

require(topGO)
require(org.Hs.eg.db)
require(KEGGprofile)
require(ggplot2)
require(biomaRt)
library(gridExtra)
library(grid)

# Initial config
config <- config::get(file = "conf/enrich_config.yml")

load(config$rdata_file)

num_clusters = 0
if (config$adt_data) {
  num_clusters = length(levels(counts_so[["rnaClusterID"]]$rnaClusterID))
} else {
  num_clusters = length(levels(counts_so@active.ident))
}

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

for (cluster_idx in 0:(config$num_clusters - 1)) {
  found_GO = TRUE
  found_KEGG = TRUE
  
  markers = markers_per_cluster[[cluster_idx + 1]]
  
  #Convert the HGNC names to Entrez IDs (eliminate non-matches where necessary)
  annots <- getBM(
    mart = mart,
    attributes = c("hgnc_symbol", "entrezgene_id"),
    filter = "hgnc_symbol",
    values = rownames(markers),
    uniqueRows = TRUE
  )
  annots <- annots[!duplicated(annots[, 1]), ]
  annots <- annots[!is.na(annots[, 2]), ]
  
  markers <- markers[which(rownames(markers) %in% annots[, 1]), ]
  annots <- annots[which(annots[, 1] %in% rownames(markers)), ]
  markers <- markers[match(annots[, 1], rownames(markers)), ]
  rownames(markers) <- annots[, 2]
  
  markers_names <- rownames(markers)
  markers <- as.numeric(markers$p_val)
  names(markers) <- markers_names
  
  #Perform gene enrichment with topGO of the target genes
  selection <- function(allScore) {
    return(allScore < 0.05)
  }
  allGO2genes <- annFUN.org(
    whichOnto = "BP",
    feasibleGenes = NULL,
    mapping = "org.Hs.eg.db",
    ID = "entrez"
  )
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = markers,
    annot = annFUN.GO2genes,
    GO2genes = allGO2genes,
    geneSel = selection,
    nodeSize = 5
  )
  
  #In order to make use of the rank information, use Kolomonogorov-Smirnov (K-S) test
  results.ks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  goEnrichment <- GenTable(GOdata,
                           KS = results.ks,
                           orderBy = "KS",
                           topNodes = 20)
  goEnrichment <- goEnrichment[goEnrichment$KS < 0.05, ]
  goEnrichment <- goEnrichment[, c("GO.ID", "Term", "KS")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep =
                               ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  
  if (dim(goEnrichment)[1] > 0) {
    p1 <- ggplot(goEnrichment, aes(x = Term, y = -log10(KS))) +
      stat_summary(
        geom = "bar",
        fill = 'royalblue',
        fun.y = mean,
        position = "dodge"
      ) +
      xlab("GO biological process") +
      ylab("Enrichment") +
      ggtitle("") +
      scale_y_continuous(breaks = round(seq(0, max(
        -log10(goEnrichment$KS)
      ), by = 2), 1)) +
      theme_bw(base_size = 24) +
      theme(
        legend.position = 'none',
        legend.background = element_rect(),
        plot.title = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          vjust = 1
        ),
        axis.text.x = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          hjust = 1.10
        ),
        axis.text.y = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          vjust = 0.5
        ),
        axis.title = element_text(size = 16, face = "bold"),
        legend.key = element_blank(),
        #removes the border
        legend.key.size = unit(1, "cm"),
        #Sets overall area/size of the legend
        legend.text = element_text(size = 12),
        #Text size
        title = element_text(size = 12)
      ) +
      guides(colour = guide_legend(override.aes = list(size = 2.5))) +
      coord_flip()
  } else {
    found_GO <- FALSE
  }
  
  #Simple pathway analysis of target genes
  markers <- markers[which(markers < 0.05)]
  bootstrap <- list()
  for (i in 1:10) {
    bootstrap[[i]] <- find_enriched_pathway(
      names(markers[sample(length(markers), length(markers) / 1.5)]),
      species = "hsa",
      returned_pvalue = 0.01,
      returned_adjpvalue = 0.05,
      returned_genenumber = 5,
      download_latest = FALSE,
      refGene = NULL
    )
    bootstrap[[i]] <- bootstrap[[i]][[1]][order(bootstrap[[i]][[1]]$pvalueAdj), ]
  }
  if (dim(do.call(rbind, bootstrap))[1] > 0) {
    keggEnrichment <- do.call(rbind, bootstrap)
    keggEnrichment <- keggEnrichment[order(keggEnrichment$pvalueAdj), ]
    keggEnrichment <- keggEnrichment[!duplicated(keggEnrichment$Pathway_Name), ]
    keggEnrichment <- keggEnrichment[keggEnrichment$pvalueAdj < 0.05, ]
    keggEnrichment$Pathway_Name <- factor(keggEnrichment$Pathway_Name,
                                          levels = rev(keggEnrichment$Pathway_Name))
    
    p2 <- ggplot(keggEnrichment, aes(x = Pathway_Name, y = -log10(pvalueAdj))) +
      stat_summary(
        geom = "bar",
        fill = 'forestgreen',
        fun.y = mean,
        position = "dodge"
      ) +
      xlab("KEGG pathway") +
      ylab("Enrichment") +
      ggtitle("") +
      scale_y_continuous(breaks = round(seq(0, max(
        -log10(keggEnrichment$pvalueAdj)
      ), by = 2), 1)) +
      theme_bw(base_size = 24) +
      theme(
        legend.position = 'none',
        legend.background = element_rect(),
        plot.title = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          vjust = 1
        ),
        axis.text.x = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          hjust = 1.10
        ),
        axis.text.y = element_text(
          angle = 0,
          size = 12,
          face = "bold",
          vjust = 0.5
        ),
        axis.title = element_text(size = 16, face = "bold"),
        legend.key = element_blank(),
        #removes the border
        legend.key.size = unit(1, "cm"),
        #Sets overall area/size of the legend
        legend.text = element_text(size = 12),
        #Text size
        title = element_text(size = 12)
      ) +
      guides(colour = guide_legend(override.aes = list(size = 2.5))) +
      coord_flip()
    
  } else {
    found_KEGG = FALSE
  }
  
  if (found_GO && found_KEGG) {
    pdf(
      paste(
        config$plots_dir,
        "cluster",
        cluster_idx,
        "_enriched_genes.pdf",
        sep = ""
      ),
      width = 18,
      height = 8.5
    )
    grid.arrange(p1, p2, ncol = 2, top = textGrob(
      paste("Cluster", cluster_idx, "enriched genes"),
      gp = gpar(fontsize = 28, fontface = "bold")
    ))
    grid.rect(gp = gpar(fill = NA))
    dev.off()
  } else if (found_GO) {
    pdf(
      paste(
        config$plots_dir,
        "cluster",
        cluster_idx,
        "_enriched_genes.pdf",
        sep = ""
      ),
      width = 9,
      height = 8.5
    )
    grid.arrange(p1, ncol = 1, top = textGrob(
      paste("Cluster", cluster_idx, "enriched genes"),
      gp = gpar(fontsize = 28, fontface = "bold")
    ))
    grid.rect(gp = gpar(fill = NA))
    dev.off()
  } else if (found_KEGG) {
    pdf(
      paste(
        config$plots_dir,
        "cluster",
        cluster_idx,
        "_enriched_genes.pdf",
        sep = ""
      ),
      width = 9,
      height = 8.5
    )
    grid.arrange(p2, ncol = 1, top = textGrob(
      paste("Cluster", cluster_idx, "enriched genes"),
      gp = gpar(fontsize = 28, fontface = "bold")
    ))
    grid.rect(gp = gpar(fill = NA))
    dev.off()
  }
}
