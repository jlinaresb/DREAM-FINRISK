setwd(here::here())
experiment_id <- "surv_xgboost"
outdir <- "src/jlb_m2/results/"

start <- Sys.time()
source("requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/jlb_m2/requirements.r")
source("src/jlb_m2/models/pipelines.r")

# Load data
train <- pseq(subset = "train")
test <- pseq(subset = "test")

# Remove samples
train <- remove_samples(train,
                        remove_nas = TRUE,
                        remove_neg = TRUE)
test <- remove_samples(test,
                       remove_nas = TRUE,
                       remove_neg = TRUE)

# Filter taxa by counts
train <- filter_taxa(train,
                     function(x) sum(x > 2) > (0.8 * length(x)), TRUE)
test <- filter_taxa(test,
                     function(x) sum(x > 2) > (0.8 * length(x)), TRUE)

# Agglomerate by Species
train <- remove_taxa(train, glomby = "Species")
train <- subset_taxa(train, Domain != "s__")
test <- remove_taxa(test, glomby = "Species")
test <- subset_taxa(test, Domain != "s__")

# Calculate richness
richness_train <- estimate_richness(
                        train,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "Simpson",
                                     "InvSimpson", "Fisher"))

richness_test <- estimate_richness(
                        test,
                        split = TRUE,
                        measures = c("Shannon", "Observed",
                                     "Chao1", "Simpson",
                                     "InvSimpson", "Fisher"))


# Generate data to train
# =======
data_train <- sample_data(train)
data_test <- sample_data(test)


data_train <- data_train[complete.cases(data_train), ]
data_test <- impute::impute.knn(t(data_test))
data_test <- as.data.frame(t(data_test$data))

# Convert variables to numeric
train_pats <- rownames(data_train)
test_pats <- rownames(data_test)
data_train <- apply(data_train, 2, function(x) as.numeric(x))
data_test <- apply(data_test, 2, function(x) as.numeric(x))

rownames(data_train) <- train_pats
rownames(data_test) <- test_pats
data_train <- as.data.frame(data_train)
data_test <- as.data.frame(data_test)

# Run training
# =====
res <- xgboost_surv_pipeline(
            data = data_train,
            dataname = experiment_id,
            time = "Event_time",
            event = "Event",
            removeConstant = TRUE,
            normalize = FALSE,
            filterFeatures = FALSE,
            inner = rsmp("holdout"),
            measure = msr("surv.cindex"),
            method_at = tnr("grid_search", resolution = 20, batch_size = 10),
            method_afs = NULL,
            term_evals = 10,
            fselector = FALSE,
            workers = 20,
            outDir = outdir,
            parallel = TRUE,
            seed = 1993
        )

preds <- res$learner$predict_newdata(data_test)

saveRDS(list(
            model = res,
            train = data_train,
            test = data_test,
            predictions = preds),
        file = "src/jlb_m2/results/survival_model.rds")


probs <- exp(-preds$crank)
time <- data_test$Event_time
event <- data_test$Event

cindex <- rcorr.cens(probs, Surv(time, event), outx = FALSE)
print(paste0("C-Index in test set model is: ", cindex[1]))
end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)