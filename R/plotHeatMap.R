plotHeatMap <- function(phyloseq, taxa_p = 0.1, seed = 0, p_adjust = "BH", group = "condition", condition_label = c("Fl", "Fl_Arg"), pvalue_cutoff = 0.05, norm = "TSS", annotation_colors = c('#FFA726','#26A69A')) {

  ps.rel = phyloseq::transform_sample_counts(phyloseq, function(x) x/sum(x)*100)
  selected_logical <<- rowSums(otu_table(ps.rel) > taxa_p) > 0
  ps.rel.filter = phyloseq::subset_taxa(phyloseq, selected_logical)

  set.seed(seed)

  ancom_output <- run_ancom(ps.rel.filter,
                            group = group,
                            taxa_rank = "all",
                            p_adjust = p_adjust,
                            pvalue_cutoff = pvalue_cutoff,
                            norm = norm)

  otu = otu_table(ancom_output)
  otu = as.data.frame(otu)

  markers = marker_table(ancom_output)
  markers = data.frame(markers)

  features <- markers$feature
  otu_selected <- otu[rownames(otu) %in% features, ]

  zscore_normalization <- function(x) {
    return((x - mean(x)) / sd(x))
  }

  markers_Fl = markers[markers$enrich_group == condition_label[1], ]
  features_Fl = markers_Fl$feature
  otu_selected_Fl = otu[rownames(otu) %in% features_Fl, ]

  markers_Fl_Arg = markers[markers$enrich_group == condition_label[2], ]
  features_Fl_Arg = markers_Fl_Arg$feature
  otu_selected_Fl_Arg = otu[rownames(otu) %in% features_Fl_Arg, ]

  otu_select_data = rbind(otu_selected_Fl, otu_selected_Fl_Arg)
  normalized_data_fa = t(apply(otu_select_data, 1, zscore_normalization))

  sample_counts <- table(sample_data(phyloseq)[[group]])
  column_annotations <- data.frame(
    Category = factor(c(rep(condition_label[1], sample_counts[condition_label[1]]), rep(condition_label[2], sample_counts[condition_label[2]])), levels = condition_label)
  )
  rownames(column_annotations) <- colnames(normalized_data_fa)

  annotation_colors <- list(
    Category = setNames(annotation_colors, condition_label)
  )

  sorted_colnames <- gtools::mixedsort(colnames(normalized_data_fa))
  normalized_data_fa <- normalized_data_fa[, sorted_colnames]

  options(repr.plot.width = 3, repr.plot.height = 6)
  pheatmap::pheatmap(normalized_data_fa,
            cluster_rows = F,
            cluster_cols = F,
            show_colnames = F,
            show_rownames = F,
            fontsize_col = 5,
            annotation_col = column_annotations,
            annotation_colors = annotation_colors,
            breaks = seq(from = -1.3, to = 1.3, length.out = 100),
            border_color = NA
  )
}
