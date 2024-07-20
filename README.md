
<p align="center">
  <img src="./Image/MicroVAR.png" alt="MicroVAR" width="200" height="200"/>
</p>

# MicroVAR

<!-- badges: start -->
<!-- badges: end -->

The MicroVAR is Microbiome Visualization and Analysis in R 

It integrates 16S rRNA data and shotgun metagenomics data using Phyloseq, generating diversity plots and taxonomy plots to visualize sample diversity and composition. Differential abundance analysis identifies significant microbial differences between groups, and microbial interaction networks are analyzed to understand the relationships between different taxa. Additionally, functional profiling with Picrust2, including KEGG annotation and pathway analysis, is performed and visualized.

- Alpha diversity measures the variety of species within a single sample. It tells us how many different species are present and how evenly they are distributed.
- Beta diversity compares the differences in species between different samples. It shows how unique or similar the species are across various environments.
- Taxonomy analysis identifies and classifies the species in a sample based on their genetic material. It helps determine which organisms are present and their abundance.
- Differential abundance analysis identifies which species or genes are present in different amounts between groups of samples. It helps in understanding how microbial communities change in response to different conditions or treatments.
- Microbial interaction networks map the relationships and interactions between different microbial species within a community. These networks reveal how species influence each other, which can help in understanding ecosystem stability and dynamics.
- Functional profiling assesses the functional capabilities of a microbial community by identifying the genes or pathways present. It provides insights into the metabolic activities and potential functions of the community, such as nutrient cycling or disease associations.

<p align="center">
<img src="./Image/Overview.png" alt="MicroVAR" width="800" height="630"/>
</p>

## Tools Comparison
<p align="center">
  <img src="./Image/tools.png" alt="MicroVAR" width="800" height="195"/>
</p>

## Learn more

* We conducted analyses using some of the raw data from the paper by [Carda-Diéguez, M., Moazzez, R. & Mira, A. Functional changes in the oral microbiome after use of fluoride and arginine containing dentifrices: a metagenomic and metatranscriptomic study. Microbiome 10, 159 (2022).](https://doi.org/10.1186/s40168-022-01338-4) 

* The raw data was downloaded from [this BioProject.](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA712952/)

* The preprocessed data used in the example can be downloaded from the link below. If you install a package, you can also load data from within the package.
  - Download the amplicon metadata from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/AmpliconMetadata.csv)
  - Download the amplicon ASV data from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/AmpliconASV.tsv)
  - Download the amplicon taxanomy data from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/AmpliconTaxonomy.tsv)
  - Download the amplicon tree data from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/AmpliconRootedTree.nwk)
  - Download the amplicon Picrust2 data from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/AmpliconPicrust.tsv)
  - Download the shotgun metadata from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/ShotgunMetadata.csv)
  - Download the shotgun biom data from [here.](https://github.com/ERASMUSlab/MicroVAR/blob/master/SampleData/ShotgunBiomData.biom)

## Installation
You can install the development version of MicroVAR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ERASMUSlab/MicroVAR")
```

## Raw Data Preprocessing
We utilized the following packages to preprocess the raw data. 
* Qiime2-2023.07
* Pircurst2 

* Trim Galore
* Bowtie2
* Samtools
* Kraken2
* Bracken  

## Begin microbiome analysis using the MicroVAR package in R
Download and load the MicroVAR package
``` r
# install.packages("devtools")
devtools::install_github("ERASMUSlab/MicroVAR")

library(MicroVAR)
```

### Creating Phyloseq Data
``` r
metadata_path <- system.file("extdata", "AmpliconMetadata.csv", package = "MicroVAR")
asv_path <- system.file("extdata", "AmpliconASV.tsv", package = "MicroVAR")
taxa_path <- system.file("extdata", "AmpliconTaxonomy.tsv", package = "MicroVAR")
tree_path <- system.file("extdata", "AmpliconRootedTree.nwk", package = "MicroVAR")

amplicon <<- createAmpPhylo(metadata_path = metadata_path,
                            asv_path = asv_path, 
                            taxa_path = taxa_path, 
                            tree_path = tree_path)
amplicon
```

```{r, message=F, warning=F, error=T, echo = TRUE}
metadata_path <- system.file("extdata", "ShotgunMetadata.csv", package = "MicroVAR")
biom_path <- system.file("extdata", "ShotgunBiomData.biom", package = "MicroVAR")

shotgun <<- createShotPhylo(biom_path = biom_path, 
                            metadata_path = metadata_path)
shotgun
```

### Creating Alpha Diversity Plot
Amplicon data
```r
shannon <- plotAlphaDiversity(phyloseq = amplicon, 
                              condition_col = "condition", 
                              condition_label = c("Fl","Fl_Arg"), 
                              colors = c('#FFA726','#26A69A'),
                              measure = "Shannon",
                              text_size = 15)
shannon

observed <- plotAlphaDiversity(phyloseq = amplicon,
                               condition_col = "condition", 
                               condition_label = c("Fl", "Fl_Arg"), 
                               colors = c('#FFA726','#26A69A'), 
                               measure ="Observed",
                               text_size = 15)
observed

simpson <- plotAlphaDiversity(phyloseq = amplicon,
                              condition_col = "condition",
                              condition_label = c("Fl", "Fl_Arg"),
                              colors = c('#FFA726','#26A69A'), 
                              measure ="Simpson",
                              text_size = 15)
simpson

chao1 <- plotAlphaDiversity(phyloseq = amplicon,  
                            condition_col = "condition", 
                            condition_label = c("Fl", "Fl_Arg"), 
                            colors = c('#FFA726','#26A69A'), 
                            measure ="Chao1",text_size = 15)

chao1 
```
Shotgun data
```r
shannon <- plotAlphaDiversity(phyloseq = shotgun, 
                              condition_col = "condition", 
                              condition_label = c("Fl","Fl_Arg"), 
                              colors = c('#FFA726','#26A69A'),
                              measure = "Shannon",
                              text_size = 15)
shannon

observed <- plotAlphaDiversity(phyloseq = shotgun,
                               condition_col = "condition", 
                               condition_label = c("Fl", "Fl_Arg"), 
                               colors = c('#FFA726','#26A69A'), 
                               measure ="Observed",
                               text_size = 15)
observed

simpson <- plotAlphaDiversity(phyloseq = shotgun,
                              condition_col = "condition",
                              condition_label = c("Fl", "Fl_Arg"),
                              colors = c('#FFA726','#26A69A'), 
                              measure ="Simpson",
                              text_size = 15)
simpson

chao1 <- plotAlphaDiversity(phyloseq = shotgun,  
                            condition_col = "condition", 
                            condition_label = c("Fl", "Fl_Arg"), 
                            colors = c('#FFA726','#26A69A'), 
                            measure ="Chao1",text_size = 15)

chao1 
```

### Creating Beta Diversity Plot
Amplicon data
```r
jaccard = plotBetaDiversity(phyloseq = amplicon,
                            condition_col = "condition", 
                            condition_label = c("Fl", "Fl_Arg"), 
                            colors = c('#FFA726','#26A69A'), 
                            measure ="jaccard",
                            text_size = 15)
jaccard

bray = plotBetaDiversity(phyloseq = amplicon, 
                         condition_col = "condition", 
                         condition_label = c("Fl", "Fl_Arg"), 
                         colors = c('#FFA726','#26A69A'), 
                         measure ="bray",
                         text_size = 15)
bray

wunifrac = plotBetaDiversity(phyloseq = amplicon,
                             condition_col = "condition", 
                             condition_label = c("Fl", "Fl_Arg"), 
                             colors = c('#FFA726','#26A69A'), 
                             measure ="wunifrac",
                             text_size = 15)
wunifrac

unifrac = plotBetaDiversity(phyloseq = amplicon, 
                            condition_col = "condition", 
                            condition_label = c("Fl", "Fl_Arg"), 
                            colors = c('#FFA726','#26A69A'), 
                            measure ="unifrac",
                            text_size = 15)
unifrac
```
Shotgun data
* The unifrac method requires a phylogenetic tree, so it cannot be used with shotgun data.
```r
jaccard = plotBetaDiversity(phyloseq = shotgun,
                            condition_col = "condition", 
                            condition_label = c("Fl", "Fl_Arg"), 
                            colors = c('#FFA726','#26A69A'), 
                            measure ="jaccard",
                            text_size = 15)
jaccard

bray = plotBetaDiversity(phyloseq = shotgun, 
                         condition_col = "condition", 
                         condition_label = c("Fl", "Fl_Arg"), 
                         colors = c('#FFA726','#26A69A'), 
                         measure ="bray",
                         text_size = 15)
bray
```
### Creating Taxonomy Plot
Amplicon
```r
kingdom_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Kingdom",
                                keep_percent = 0,
                                text_size = 10)
kingdom_taxa

phylum_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Phylum",
                                keep_percent = 1,
                                text_size = 10)
phylum_taxa

class_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Class",
                                keep_percent = 1,
                                text_size = 10)
class_taxa

order_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Order",
                                keep_percent = 1,
                                text_size = 10)
order_taxa

family_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Family",
                                keep_percent = 1,
                                text_size = 7)
family_taxa

genus_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                               condition_label = c("Fl","Fl_Arg"),
                               taxa_level = "Genus",
                               keep_percent = 1,
                               text_size = 7)
genus_taxa

species_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                 condition_label = c("Fl","Fl_Arg"),
                                 taxa_level = "Species",
                                 keep_percent = 1,
                                 text_size = 10)
species_taxa
```
Shotgun
```r
kingdom_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Kingdom",
                                keep_percent = 0,
                                text_size = 10)
kingdom_taxa

phylum_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Phylum",
                                keep_percent = 1,
                                text_size = 10)
phylum_taxa

class_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Class",
                                keep_percent = 1,
                                text_size = 10)
class_taxa

order_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Order",
                                keep_percent = 1,
                                text_size = 10)
order_taxa

family_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                                condition_label = c("Fl","Fl_Arg"),
                                taxa_level = "Family",
                                keep_percent = 1,
                                text_size = 10)
family_taxa

genus_taxa <- plotTaxonomyBars(phyloseq = shotgun,
                               condition_label = c("Fl","Fl_Arg"),
                               taxa_level = "Genus",
                               keep_percent = 1,
                               text_size = 7)
genus_taxa

species_taxa <- plotTaxonomyBars(phyloseq = amplicon,
                                 condition_label = c("Fl","Fl_Arg"),
                                 taxa_level = "Species",
                                 keep_percent = 1,
                                 text_size = 6)
species_taxa
```
### Differential abundance plot
Heat tree
```r
plotHeatTree(phyloseq = amplicon,
             taxa_p = 0.1,
             seed = 0,
             p_adjust = "BH",
             group ="condition",
             condition_label = c("Fl", "Fl_Arg"),
             pvalue_cutoff = 0.05, 
             norm = "TSS",
             clade_label_level = 7,
             lade_label_font_size = 10,
             colors = c('#FFA726','#26A69A'),
             text_size = 10)
```
Heat map 
```r
plotHeatMap(phyloseq = amplicon,
            taxa_p = 0.1,
            seed = 0,
            p_adjust = "BH",
            group = "condition",
            condition_label = c("Fl", "Fl_Arg"), 
            pvalue_cutoff = 0.05,
            norm = "TSS",
            annotation_colors = c('#FFA726', '#26A69A'))
```



