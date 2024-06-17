createAmpPhylo <- function(asv_path, taxa_path, metadata_path, tree_path) {

  options(warn=-1)

  # taxanomy 결과 정리 함수
  clean_taxonomy <- function(taxa_path) {

    # 데이터 읽기
    taxonomy = read.table(taxa_path, sep = "\t", header = TRUE, row.names = 1)

    # Taxon 분리
    tax = taxonomy %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")

    tax_clean = data.frame(row.names = row.names(tax),
                           Kingdom = str_replace(tax[,1], "d__",""),
                           Phylum = str_replace(tax[,2], "p__",""),
                           Class = str_replace(tax[,3], "c__",""),
                           Order = str_replace(tax[,4], "o__",""),
                           Family = str_replace(tax[,5], "f__",""),
                           Genus = str_replace(tax[,6], "g__",""),
                           Species = str_replace(tax[,7], "s__",""),
                           stringsAsFactors = FALSE)

    tax_clean[is.na(tax_clean)] = ""
    tax_clean[tax_clean == "__"] = ""

    for (i in 1:nrow(tax_clean)){
      if (tax_clean[i,7] != ""){
        tax_clean$Species[i] = paste("", tax_clean$Species[i], sep = "")
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
        tax_clean$Species[i] = paste("Unclassified ",tax_clean$Genus[i], sep = " ")
      }
    }
    return(tax_clean)
  }

  # taxonomy 정리 함수 활용
  asv = read.table(asv_path,
                   sep = "\t",header=TRUE,row.names = 1)
  clean_taxa = clean_taxonomy(taxa_path)

  # 샘플 데이터
  metadata = read.csv(metadata_path, row.names = 1)

  # phyloseq 객체 생성
  ASV = otu_table(as.matrix(asv), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(clean_taxa))
  SAMPLE = sample_data(metadata)
  TREE = read_tree(tree_path)
  amplicon <<- phyloseq::phyloseq(ASV, TAX, SAMPLE, TREE)

  return(amplicon)
}
