plotTaxonomyBars <- function(phyloseq, condition_label, taxa_level, keep_percent,text_size = 20) {

  ps.rel = transform_sample_counts(phyloseq, function(x) x/sum(x)*100)

  glom = tax_glom(ps.rel, taxrank = taxa_level , NArm = FALSE)
  ps.melt = psmelt(glom)

  ps.melt[[taxa_level]] = as.character(ps.melt[[taxa_level]])

  ps.melt = ps.melt %>%
    group_by(condition, !!sym(taxa_level)) %>%
    mutate(median=median(Abundance))

  keep_label = paste0("< ", keep_percent, "%")

  keep = unique(ps.melt[[taxa_level]][ps.melt$median > keep_percent])
  ps.melt[[taxa_level]][!(ps.melt[[taxa_level]] %in% keep)] = keep_label


  ps.melt_sum <- ps.melt %>%
    group_by(Sample, condition, !!sym(taxa_level)) %>%
    summarise(Abundance = sum(Abundance))

  means = aggregate(Abundance ~ condition + ps.melt_sum[, 3, drop = TRUE], data = ps.melt_sum, FUN = mean)
  means$condition = factor(means$condition, levels = condition_label)



  p = ggplot(means, aes(x = condition, y = Abundance, fill = means[,2])) +
    geom_bar(stat = "identity", position = "stack",color='black') +
    labs(x="", y="Relative Abundance(%)",fill=taxa_level) +
    facet_wrap(~condition, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(strip.background = element_blank(),
          axis.text.x.bottom = element_text(size = text_size, color = "black", face = "bold"),
          axis.text.y = element_text(size = text_size, color = "black", face = "bold"),
          axis.text.x = element_text(size = text_size, color = "black", face = "bold"),
          axis.title = element_text(size = text_size, color = "black", face = "bold"),
          legend.text = element_text(size = text_size, color = "black", face = "bold"),
          legend.title = element_text(size = text_size, color = "black", face = "bold"),
          strip.text = element_blank())




  return(p)
}
