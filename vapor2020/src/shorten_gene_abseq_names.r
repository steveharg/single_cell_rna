shorten_gene_abseq_names <- function(counts_df) {
  orig_gene_abseq_names <- rownames(counts_df)
  gene_abseq_name_splits <- strsplit(rownames(counts_df), "\\.")
  gene_abseq_name_split_lengths <- sapply(X = gene_abseq_name_splits, FUN = length)
  short_gene_abseq_names <- vector("list", length(gene_abseq_name_splits))
  
  for (i in 1:length(gene_abseq_name_splits)) {
    is_protein <- gene_abseq_name_splits[[i]][length(gene_abseq_name_splits[[i]])] == "pAbO"
    
    if (is_protein) {
      if (gene_abseq_name_split_lengths[[i]] == 4) {
        short_gene_abseq_names[[i]] <- paste0(
          gene_abseq_name_splits[[i]][1],
          ".",
          gene_abseq_name_splits[[i]][2],
          ".",
          gene_abseq_name_splits[[i]][4]
        )
      }
      else {
        short_gene_abseq_names[[i]] <- paste0(
          gene_abseq_name_splits[[i]][1],
          ".",
          gene_abseq_name_splits[[i]][2],
          ".",
          gene_abseq_name_splits[[i]][3],
          ".",
          gene_abseq_name_splits[[i]][5]
        )
      }
    }
    else {
      if (gene_abseq_name_split_lengths[[i]] == 4) {
        short_gene_abseq_names[[i]] <- gene_abseq_name_splits[[i]][1]
      }
      else {
        short_gene_abseq_names[[i]] <- paste0(gene_abseq_name_splits[[i]][1],
                                              ".",
                                              gene_abseq_name_splits[[i]][2])
      }
    }
  }
  
  return(short_gene_abseq_names)
}