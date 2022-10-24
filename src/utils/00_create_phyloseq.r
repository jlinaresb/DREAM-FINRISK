# TODO:
# - Remove non-targeting patients (see discussion forum)
# - Add guanrank variable to stratify patients
# - Remove NZV variables
# - Calculate relative abundance
# - Feature selection for each taxon level (PCA, FCBF, filter, wrapper)
# - Trainning models:
#     Option 1): Survival + treatment (like anthra project)
#     Option 2): Binary classification after guanrank stratification
#     Option 3): Survival model 


# Create phyloseq object
# ===
require(phyloseq)
setwd("~/git/DREAM-FINRISK/extdata/DreamHF/train/")

# Clinical
pheno = read.csv("pheno_training.csv", header = T, row.names = 1)

# Tax table
taxtable = read.csv("taxtable.csv", header = T)
rownames(taxtable) = paste("taxid", rownames(taxtable), sep = "_")
taxtable = as.matrix(taxtable)

# OTUs
counts = read.csv("readcounts_training.csv", header = T, row.names = 1)
rownames(counts) = rownames(taxtable)
counts = as.matrix(counts)

# Create phyloseq
OTU = otu_table(counts, taxa_are_rows = TRUE)
TAX = tax_table(taxtable)
SAMPLE = sample_data(pheno)

pseq = phyloseq(OTU, TAX, SAMPLE)
