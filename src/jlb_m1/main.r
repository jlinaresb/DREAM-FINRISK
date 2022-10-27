# TODO:
# - Glomerate by species/genus?
# - Feature reduction for each taxon level (PCA, FCBF)
# - Return scores.csv
# - Do not remove any sample in test!


setwd(here::here())
experiment_id <- "pca"
outdir <- "src/jlb_m1/results/"

start <- Sys.time()
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/jlb_m1/fit_model.r")
source("src/jlb_m1/predRes_helper.r")

# Import data
# ======
train <- pseq(subset = "train")
test <- pseq(subset = "test")


# Preprocess data
# ======
# Remove samples
train <- remove_samples(train,
                        remove_nas = TRUE,
                        remove_neg = FALSE)
test <- remove_samples(test,
                       remove_nas = TRUE,
                       remove_neg = FALSE)

# Relabel patients with PrevalentHFAIL
sample_data(train)[sample_data(train)$PrevalentHFAIL == 1 &
                   sample_data(train)$Event == 0, ]$Event_time <-
                   max(sample_data(train)$Event_time)
sample_data(test)[sample_data(test)$PrevalentHFAIL == 1 &
                   sample_data(test)$Event == 0, ]$Event_time <-
                   max(sample_data(test)$Event_time)

# Calculate richness
richness_train <- estimate_richness(
                        train,
                        split = TRUE,
                        measures = c("Shannon"))
sample_data(train)$shannon_index <- richness_train$Shannon

richness_test <- estimate_richness(
                        test,
                        split = TRUE,
                        measures = c("Shannon"))
sample_data(test)$shannon_index <- richness_test$Shannon

# Remove and normalization taxa in train
train <- remove_taxa(train)
train <- norm(train)
train <- subset_taxa(train, Species != "s__")

test <- norm(test, filter = FALSE)
taxids <- taxa_names(train)
test <- prune_taxa(taxids, test)


# PCA
# ======
otu_train <- t(otu_table(train)@.Data)
otu_test <- t(otu_table(test)@.Data)
pca <- prcomp(otu_train, scale = TRUE)
otu_train_pca <- pca$x[, 1:10]
otu_test_pca <- predict(pca, newdata = otu_test)[, 1:10]


# Generate data to train
# =======
pheno_train <- sample_data(train)
pheno_test <- sample_data(test)
data_train <- cbind.data.frame(pheno_train, otu_train_pca)
data_test <- cbind.data.frame(pheno_test, otu_test_pca)

data_train <- data_train[complete.cases(data_train), ]
data_test <- data_test[complete.cases(data_test), ]


# Fit model
# ========
methods <- c("alassoR", "alassoU",
             "enet", "lasso",
             "lasso-pcvl", "lasso-RIC",
             "PCAlasso", "ridgelasso")
models <- fit_biospear(data = data_train,
                       biomarkers = paste0("PC", 1:10),
                       surv = c("Event_time", "Event"),
                       cvrts = c("Age", "BodyMassIndex", "Smoking",
                                 "PrevalentDiabetes", "PrevalentCHD",
                                 "SystolicBP", "NonHDLcholesterol",
                                 "Sex", "shannon_index"),
                       treatment = "BPTreatment",
                       methods = methods)

environment(predRes2) <- asNamespace('biospear')
assignInNamespace("predRes", predRes2, ns = "biospear")
prediction <- predRes2(res = models,
                    method = methods,
                    traindata = data_train,
                    newdata = data_test,
                    int.cv = FALSE,
                    int.cv.nfold = 5,
                    time = seq(1, 16, 1),
                    trace = TRUE,
                    ncores = 8)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)


# Save data and models
# =======
saveRDS(list(
            train = data_train,
            test = data_test,
            models = models,
            prediction = prediction),
            file = paste0(outdir, "model_", experiment_id, ".rds"))
