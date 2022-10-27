# Function to import train or test data in phyloseq object
# ========
# Arguments:
#   - subset: character indicating "train" or "test" subset.
# Return:
#   phyloseq object with pheno, raw counts and taxtable


pseq <- function(subset) {
  folder_files <- list.files(path = subset)
  # Clinical
  pheno_file <- folder_files[grepl("pheno", folder_files)]
  pheno <- read.csv(paste0(subset, "/", pheno_file),
                      header = TRUE, row.names = 1)
  # Tax table
  taxtable <- read.csv(paste0(subset, "/taxtable.csv"),
                         header = TRUE)
  rownames(taxtable) <- paste("taxid", rownames(taxtable), sep = "_")
  taxtable <- as.matrix(taxtable)

  # OTUs
  counts_file <- folder_files[grepl("counts",folder_files)]
  counts <- read.csv(paste0(subset, "/",counts_file),
                      header = TRUE, row.names = 1)
  rownames(counts) <- rownames(taxtable)
  counts <- as.matrix(counts)

  # Create phyloseq
  otu <- otu_table(counts, taxa_are_rows = TRUE)
  tax <- tax_table(taxtable)
  sample <- sample_data(pheno)

  pseq <- phyloseq(otu, tax, sample)

  return(pseq)
}