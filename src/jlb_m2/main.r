setwd(here::here())
experiment_id <- "diff_rf"
outdir <- "src/jlb_m2/results/"

star <- Sys.time()
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/utils/guanrank.r")
source("src/utils/dea.r")
source("src/jlb_m2/requirements.r")
source("src/jlb_m2/models/pipelines.r")


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

# Calculate guanrank
gr <- data.frame(
    time = data_train$Event_time,
    status = data_train$Event,
    row.names = rownames(data_train))
gr <- as.data.frame(guanrank(gr))

low_risk <- gr[which(gr$rank < 0.45), ]
high_risk <- gr[which(gr$rank > 0.7), ]

low_risk$target <- "low_risk"
high_risk$target <- "high_risk"
data <- rbind.data.frame(high_risk, low_risk)
data_train <- data_train[match(rownames(data), rownames(data_train)), ]

to_train <- data.frame(
    subset(data_train, select = -c(Event, Event_time)),
    target = data$target,
    row.names = rownames(data_train)
)

# Run training
# =====
rf_pipeline(
    data = to_train,
    dataname = experiment_id,
    target = "target",
    positive = "high_risk",
    removeConstant = TRUE,
    normalize = TRUE,
    filterFeatures = FALSE,
    inner = rsmp("holdout", ratio = 0.7),
    outer = rsmp("cv", folds = 10),
    measure = msr("classif.auc"),
    method_at = tnr("grid_search", resolution = 30, batch_size = 10),
    method_afs = NULL,
    term_evals = NULL,
    fselector = FALSE,
    workers = 20,
    outDir = outdir,
    parallel = TRUE,
    seed = 1993
)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)