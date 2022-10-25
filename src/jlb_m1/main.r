# TODO:
# - Glomerate by species/genus?
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
saveRDS(list(
            train = data,
            fit = models),
            file = "src/jlb_m1/results/models.rds")


# Validation
# ========
test <- pseq(subset = "test")

test <- remove_taxa(test)
test <- norm(test)

taxids <- taxa_names(train)
test <- prune_taxa(taxids, test)

otu_test <- as.data.frame(t(otu_table(test)@.Data))
pheno_test <- sample_data(test)
data_test <- cbind.data.frame(pheno_test, otu_test)

predRes(res = models,
        method = methods,
        traindata = data,
        newdata = test,
        int.cv = TRUE,
        int.cv.nfold = 5,
        time = seq(2, 15, 1),
        trace = TRUE,
        ncores = 5
        )