createShotPhylo <- function(biom_path, metadata_path) {

  options(warn=-1)

  # taxanomy 결과 정리 함수
  clean_taxonomy <- function(ps) {

    # 데이터 읽기
    taxonomy  = tax_table(shotgun)

    tax_clean = data.frame(row.names = row.names(taxonomy),
                           Kingdom = str_replace(taxonomy[,1], "k__",""),
                           Phylum = str_replace(taxonomy[,2], "p__",""),
                           Class = str_replace(taxonomy[,3], "c__",""),
                           Order = str_replace(taxonomy[,4], "o__",""),
                           Family = str_replace(taxonomy[,5], "f__",""),
                           Genus = str_replace(taxonomy[,6], "g__",""),
                           Species = str_replace(taxonomy[,7], "s__",""),
                           stringsAsFactors = FALSE)

    tax_clean[is.na(tax_clean)] = ""
    tax_clean[tax_clean == "__"] = ""

    for (i in 1:nrow(tax_clean)){
      if (tax_clean[i,7] != ""){
        tax_clean$Species[i] = paste(tax_clean$Genus[i], tax_clean$Species[i], sep = " ")
      } else if (tax_clean[i,2] == ""){
        kingdom = paste("Unclassified", tax_clean[i,1], sep = " ")
        tax_clean[i, 2:7] = kingdom
      } else if (tax_clean[i,3] == ""){
        phylum = paste("Unclassified", tax_clean[i,2], sep = " ")
        tax_clean[i, 3:7] = phylum
      } else if (tax_clean[i,4] == ""){
        class = paste("Unclassified", tax_clean[i,3], sep = " ")
        tax_clean[i, 4:7] = class
      } else if (tax_clean[i,5] == ""){
        order = paste("Unclassified", tax_clean[i,4], sep = " ")
        tax_clean[i, 5:7] = order
      } else if (tax_clean[i,6] == ""){
        family = paste("Unclassified", tax_clean[i,5], sep = " ")
        tax_clean[i, 6:7] = family
      } else if (tax_clean[i,7] == ""){
        tax_clean$Species[i] = paste("Unclassified",tax_clean$Genus[i], sep = " ")
      }
    }
    return(tax_clean)
  }

  shotgun <<- import_biom(biom_path)

  # taxonomy 정리 함수 활용
  otu = otu_table(shotgun)
  colnames(otu) = gsub("_bracken_species", "", colnames(otu))
  clean_taxa = clean_taxonomy(shotgun)

  # 샘플 데이터
  metadata = read.csv(metadata_path, row.names = 1)

  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(clean_taxa))
  SAMPLE = sample_data(metadata)

  shotgun <<- phyloseq(OTU, TAX, SAMPLE)

  return(shotgun)
}
