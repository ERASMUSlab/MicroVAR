plotHeatMap <- function(phyloseq, taxa_p, seed, p_adjust = "BH", group, condition_label, pvalue_cutoff = 0.05, norm = "TSS",
                        colors = c('#F06292','#4DD0E1'),cluster_rows = F, cluster_cols = T, show_colnames = F, show_rownames = F ) {

  ps.rel <- transform_sample_counts(phyloseq, function(x) x/sum(x)*100)
  selected_logical <- rowSums(otu_table(ps.rel) > taxa_p) > 0
  ps.rel.filter <- subset_taxa(phyloseq, selected_logical)

  set.seed(seed)

  ancom_output <- run_ancom(ps.rel.filter,
                            group = group,
                            taxa_rank = "all",
                            p_adjust = p_adjust,
                            pvalue_cutoff = pvalue_cutoff,
                            norm = norm)

  otu <- otu_table(ancom_output)
  otu <- data.frame(otu)

  markers <- marker_table(ancom_output)
  markers <- data.frame(markers)

  features <- markers$feature
  otu_selected <- otu[rownames(otu) %in% features, ]

  zscore_normalization <- function(x) {
    return((x - mean(x)) / sd(x))
  }

  markers_group1 <- markers[markers$enrich_group == condition_label[1], ]
  features_group1 <- markers_group1$feature
  otu_selected_group1 <- otu[rownames(otu) %in% features_group1, ]

  markers_group2 <- markers[markers$enrich_group == condition_label[2], ]
  features_group2 <- markers_group2$feature
  otu_selected_group2 <- otu[rownames(otu) %in% features_group2, ]

  otu_select_data <- rbind(otu_selected_group1, otu_selected_group2)
  normalized_data_fa <- t(apply(otu_select_data, 1, zscore_normalization))

  column_annotations <- data.frame(
    Category = c(rep(condition_label[1], ncol(otu_selected_group1)), rep(condition_label[2], ncol(otu_selected_group2)))
  )
  rownames(column_annotations) <- colnames(normalized_data_fa)

  annotation_colors <- list(
    Category = setNames(colors, condition_label)
  )


  options(repr.plot.width = 3, repr.plot.height = 6)
  pheatmap(
    normalized_data_fa,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_colnames = show_colnames,
    show_rownames = show_rownames,
    annotation_col = column_annotations,
    annotation_colors = annotation_colors,
    breaks = seq(from = -2, to = 2, length.out = 100)
  )
}
