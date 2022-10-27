setwd(here::here())
experiment_id <- "diff_exp"
outdir <- "src/jlb_m1/results/"

start <- Sys.time()
source("requirements.r")
source("src/utils/dea.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/jlb_m1/fit_model.r")
source("src/jlb_m1/predRes_helper.r")

# Load data
train <- pseq(subset = "train")
test <- pseq(subset = "test")


# Differential Expression
de <- subset_samples(train, PrevalentHFAIL == 0 &
                               Smoking == 0 &
                               BPTreatment == 0 &
                               PrevalentDiabetes == 0 &
                               PrevalentCHD == 0 &
                               BodyMassIndex < 30 &
                               BodyMassIndex > 18)

de <- remove_taxa(de, glomby = "Genus")
de <- subset_taxa(de, Genus != "g__")

sample_data(de)$Event <- as.factor(sample_data(de)$Event)
sig_tab <- phyloseq_dea(pseq = de,
                    test = "Wald",
                    fit_type = "parametric",
                    alpha = 0.05)
taxids <- rownames(sig_tab)


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


# Prune taxa according differential expression
train <- prune_taxa(taxids, train)
test <- prune_taxa(taxids, test)

# Remove samples
train <- remove_samples(train,
                        remove_nas = TRUE,
                        remove_neg = FALSE)
test <- remove_samples(test,
                       remove_nas = TRUE,
                       remove_neg = FALSE)

# Generate data to train
# =======
otu_train <- apply(t(otu_table(train)@.Data), 2, function(x) log2(x + 1))
otu_test <- apply(t(otu_table(test)@.Data), 2, function(x) log2(x + 1))

pheno_train <- sample_data(train)
pheno_test <- sample_data(test)

data_train <- cbind.data.frame(pheno_train, otu_train)
data_test <- cbind.data.frame(pheno_test, otu_test)

data_train <- data_train[complete.cases(data_train), ]
data_test <- impute::impute.knn(t(data_test))
data_test <- as.data.frame(t(data_test$data))

# Relabel patients with PrevalentHFAIL
data_train$Event_time[data_train$Event_time < 0] <- 15
data_test$Event_time[data_test$Event_time < 0] <- 15


# Fit model
# ========
methods <- c("alassoR", "alassoU",
             "enet", "lasso",
             "lasso-pcvl", "lasso-RIC",
             "PCAlasso", "ridgelasso")
models <- fit_biospear(data = data_train,
                       biomarkers = grep("taxid", colnames(data_train)),
                       surv = c("Event_time", "Event"),
                       cvrts = c("Age", "BodyMassIndex", "Smoking",
                                 "PrevalentDiabetes", "PrevalentCHD",
                                 "SystolicBP", "NonHDLcholesterol",
                                 "Sex", "shannon_index", "BPTreatment"),
                       methods = methods)


environment(predRes2) <- asNamespace("biospear")
assignInNamespace("predRes", predRes2, ns = "biospear")
prediction <- predRes2(res = models,
                    method = methods,
                    traindata = data_train,
                    newdata = data_test,
                    int.cv = TRUE,
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