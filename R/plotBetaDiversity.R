plotBetaDiversity <- function(phyloseq, condition_col, condition_label, colors = c('#4DD0E1','#F06292'), measure,text_size = 20) {

  sample_data_df = phyloseq::sample_data(phyloseq)
  sample_data_df[[condition_col]] = factor(sample_data_df[[condition_col]], levels = condition_label)
  phyloseq::sample_data(phyloseq) = sample_data_df

  dist = phyloseq::distance(phyloseq, method=measure)
  ordination = ordinate(phyloseq, method="PCoA", distance=dist)

  Title = ""
  if (measure=="bray"){
    Title = "Bray-Curtis"
  } else if (measure=="jaccard"){
    Title = "Jaccard Distance"
  } else if (measure=="wunifrac"){
    Title = "Weighted Unifrac"
  } else if (measure=="unifrac"){
    Title = "Unifrac"
  }

  p = plot_ordination(phyloseq, ordination, color="condition") +
    theme_classic() +
    geom_point(size = 4) +
    scale_color_manual(values=colors, breaks = condition_label) +
    stat_ellipse(geom = "polygon", aes(color = condition), alpha = 0.03) +
    labs(title=Title) +
    theme(strip.background = element_blank(),
          axis.text.y = element_text(size = text_size, color = "black", face = "bold"),
          axis.text.x = element_text(size = text_size, color = "black", face = "bold"),
          axis.title = element_text(size = text_size, color = "black", face = "bold"),
          legend.text = element_text(size = text_size, color = "black", face = "bold"),
          legend.title = element_text(size = text_size, color = "black", face = "bold"),
          plot.title = element_text(size = text_size * 1.5, color = "black", face = "bold", hjust = 0.5))

  return(p)
}
