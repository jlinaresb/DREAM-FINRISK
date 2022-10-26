# TODO:
# - Glomerate by species/genus?
# - Calculate richness scores
# - GSEA/GSVA for microbes sets
# - Feature selection for each taxon level (PCA, FCBF)
# - Trainning models:
#     Option 1): Survival + treatment (like anthra project)
# - Process test data

start <- Sys.time()
setwd(here::here())
source("requirements.r")
source("src/utils/importPseq.r")
source("src/jlb_m1/preprocessing.r")
source("src/jlb_m1/fit_model.r")
source("src/jlb_m1/predRes_modified.r")

# Import data
# ======
train <- pseq(subset = "train")

# Preprocess data
# ======
# Remove samples
train <- remove_samples(train)
# Calculate richness
richness <- estimate_richness(train,
                        split = TRUE,
                        measures = c("Shannon", "Chao1"))

# Remove taxa
train <- remove_taxa(train)
# Normalization
train <- norm(train)
# Remove unknown species
train <- subset_taxa(train, Species != "s__")



# Fit model
# ======
otu <- as.data.frame(t(otu_table(train)@.Data))
pheno <- sample_data(train)
data <- cbind.data.frame(pheno, otu)

methods <- c("alassoR", "alassoU",
             "enet",  "lasso",
             "lasso-pcvl", "lasso-RIC",
             "PCAlasso", "ridgelasso")

models <- fit_biospear(data = data,
                       biomarkers = c(13:205),
                       surv = c("Event_time", "Event"),
                       cvrts = c(1:3, 5:6, 10:12),
                       treatment = "BPTreatment",
                       methods = methods)

end <- Sys.time()
time <- difftime(end, start, units = "hours")
print(time)

# Save train data and models
saveRDS(list(
            train = data,
            fit = models),
            file = "src/jlb_m1/results/models.rds")