setwd(here::here())
outdir <- "src/jlb_m2/results/"

source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/utils/dea.r")

# Load train
train <- pseq(subset = "train")

# Agglomerate
train <- tax_glom(train, "Genus")
# train <- subset_taxa(train, Species != "s__")

# Filter taxa
# train <- filter_taxa(train,
#                         function(x) sum(x > 5) > (0.3 * length(x)), TRUE)

otus <- otu_table(train)@.Data
pheno <- sample_data(train)
pheno <- pheno[, c("Age", "BodyMassIndex", "SystolicBP", "NonHDLcholesterol")]

# Get correlations
corr <- cor(pheno, t(otus), use = "everything", method = "pearson")
corr[1:4, 1:4]
dim(corr)
hist(as.vector(corr[3, ]))
