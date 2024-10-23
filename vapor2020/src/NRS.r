#Input: a counts.csv file and a cells to clusters mapping csv file
#Output: .txt files of NRS scores for each marker, per cluster and for the entire dataset

library(Seurat)

# Initial config
config_tmp <- config::get(file = "conf/nrs_config.yml")

# Function for calculating non-redundancy score, from
# Nowicka, Malgorzata, et al. "CyTOF workflow: differential discovery in high-throughput high-dimensional cytometry datasets."
# F1000Research 6 (2017).
NRS <- function(x, ncomp = 2) {
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  # proportionvar <- ((pr$sdev ^ 2) / (sum(pr$sdev ^ 2))) * 100
  score <- rowSums(outer(rep(1, ncol(x)), pr$sdev[1:ncomp] ^ 2) * abs(pr$rotation[, 1:ncomp]))
  return(score)
}

load(config_tmp$rdata_file)

config = config_tmp
rm(config_tmp)

counts = t(eval(parse(text = config$counts_object)))

# Indentify which cells belong to which cluster, and the counts for those cells
counts_cluster = list(num_clusters)
for (i in 1:num_clusters) {
  counts_cluster[[i]] <- t(data.frame(GetAssayData(counts_so[, cluster_cells[[i]]], slot =
                                                     "data")))
}

# Calculate NRS for the entire dataset
nrs_all = lapply(list(counts), NRS)

# Calculate NRS for the each cluster
nrs_per_cluster = lapply(counts_cluster, NRS)

# Write sorted results to files
write.table(sort(nrs_all[[1]], decreasing = T),
            file = 'data/nrs_all.txt',
            col.names = FALSE)

for (i in 1:num_clusters) {
  write.table(
    sort(nrs_per_cluster[[i]], decreasing = T),
    file = paste('data/nrs_cluster_', (i - 1), '.txt', sep = ''),
    col.names = FALSE
  )
}
