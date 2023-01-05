# Function to import train or test data in phyloseq object
# ========
# Arguments:
#   - inputdir: character indicating "train" or "test" subset.
# Return:
#   phyloseq object with pheno, raw counts and taxtable


# Utils functions
pseq <- function(inputdir, subset) {
  dir <- file.path(inputdir, subset)
  folder_files <- list.files(path = dir)
  # Clinical
  pheno_file <- folder_files[grepl("pheno", folder_files)]
  pheno <- read.csv(file.path(dir, "/", pheno_file), row.names = 1)
  # Tax table
  taxtable <- read.csv(file.path(dir, "/taxtable.csv"))
  rownames(taxtable) <- paste("taxid", rownames(taxtable), sep = "_")
  taxtable <- as.matrix(taxtable)

  # OTUs
  counts_file <- folder_files[grepl("counts", folder_files)]
  counts <- read.csv(file.path(dir, counts_file), row.names = 1)
  rownames(counts) <- rownames(taxtable)
  counts <- as.matrix(counts)

  # Create phyloseq
  otu <- otu_table(counts, taxa_are_rows = TRUE)
  tax <- tax_table(taxtable)
  sample <- sample_data(pheno)

  pseq <- phyloseq(otu, tax, sample)

  return(pseq)
}