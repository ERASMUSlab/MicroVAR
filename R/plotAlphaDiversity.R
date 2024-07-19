plotAlphaDiversity <- function(phyloseq, condition_col, condition_label, colors = c('#4DD0E1','#F06292'), measure,text_size = 20) {

  sample_data_df = phyloseq::sample_data(phyloseq)
  sample_data_df[[condition_col]] = factor(sample_data_df[[condition_col]], levels = condition_label)
  phyloseq::sample_data(phyloseq) = sample_data_df

  Title = ""
  if (measure=="Shannon"){
    Title = "Shannon Index"
  } else if (measure=="Observed"){
    Title = "Observed Feature"
  } else if (measure=="Simpson"){
    Title = "Simpson Index"
  } else if (measure=="Chao1"){
    Title = "Chao1"
  }

  p = phyloseq::plot_richness(phyloseq, x="condition",color = "condition", measures=c(measure)) +
      theme_classic() +
      facet_null()+
      geom_violin(adjust=0.7, scale='width', fill = "grey", alpha=0.1) +
      geom_boxplot(width=0.3) +
      geom_jitter(color="black", size=1, alpha=0.8) +
      geom_point(size = 4) +
      scale_color_manual(values=colors, breaks = condition_label) +
      labs(x = "Condition", y = 'Index value', title=Title) +
      theme(strip.background = element_blank(),
            axis.text.x = element_text(size = text_size, color = "black", face = "bold"),
            axis.text.y = element_text(size = text_size, color = "black", face = "bold"),
            axis.title = element_text(size = text_size, color = "black", face = "bold"),
            legend.text = element_text(size = text_size, color = "black", face = "bold"),
            legend.title = element_text(size = text_size, color = "black", face = "bold"),
            plot.title = element_text(size = text_size * 1.5, color = "black", face = "bold", hjust = 0.5))

  return(p)
}
