plotHeatTree <- function(phyloseq, taxa_p, seed, p_adjust = "BH", group, condition_label ,pvalue_cutoff = 0.05, norm ="TSS",
                           colors = c('#F06292','#4DD0E1'),clade_label_level = 7,lade_label_font_size = 6, text_size = 15) {

  ps.rel = phyloseq::transform_sample_counts(phyloseq, function(x) x/sum(x)*100);
  selected_logical <<- rowSums(otu_table(ps.rel) > taxa_p) > 0;
  ps.rel.filter = phyloseq::subset_taxa(phyloseq, selected_logical);

  set.seed(seed)

  ancom_output <- run_ancom(ps.rel.filter,
                            group = group,
                            taxa_rank = "all",
                            p_adjust = p_adjust,
                            pvalue_cutoff = 0.05,
                            norm = norm)


  options(repr.plot.width = 30, repr.plot.height = 20)
  p = plot_cladogram(ancom_output, color = colors, clade_label_level = clade_label_level, clade_label_font_size=lade_label_font_size)
  p = p + theme(strip.background = element_blank(),
                legend.text = element_text(size = text_size, color = "black", face = "bold"),
                legend.title = element_text(size = text_size, color = "black", face = "bold"),
                legend.position = "bottom")

  return(p)
}
