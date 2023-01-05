# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_2"
inputdir <- args[1]
outputdir <- file.path(inputdir, Team_Name_Submission_Number, "output")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)

start <- Sys.time()
source("/mnt/utils/requirements.r")
source("/mnt/utils/importPseq.r")
source("/mnt/utils/prepro_functions.r")
source("/mnt/utils/get_scores.r")
source("/mnt/utils/co-abundances.r")

source("/mnt//models/tuning_utils.r")
source("/mnt/models/build_learners.r")
source("/mnt/models/pipeline_utils.r")
source("/mnt/models/pipelines.r")

# Preprocess
# =======
source("/mnt/utils/preprocessing.r")

# Check data
# ======
stopifnot(ncol(x_train) == ncol(x_test))
stopifnot(colnames(x_train) == colnames(x_test))
stopifnot(nrow(x_train) == nrow(y_train))
stopifnot(nrow(x_test) == nrow(y_test))

# Creating train and test data
print("Creating train and test data ...")
train <- cbind.data.frame(x_train, y_train)
test  <- cbind.data.frame(x_test, y_test)

# What do we do with ...
# ======
# Negatives survival values in train and test?
train <- train[-which(train$Event_time < 0), ]
test$Event_time[which(test$Event_time < 0)] <- 15

# NA's values in train and test?
train <- train[complete.cases(train), ]
test <- missRanger(test, pmm.k = 10, seed = 153)

# Removing PrevalentHFAIL (only 0)
train <- train[, -grep("PrevalentHFAIL", colnames(train))]
test <- test[, -grep("PrevalentHFAIL", colnames(test))]


# Run training
# =====
res <- xgboost_surv_pipeline(
            data = train,
            dataname = "xgboost_survival",
            time = "Event_time",
            event = "Event",
            removeConstant = TRUE,
            normalize = FALSE,
            filterFeatures = FALSE,
            inner = rsmp("cv", folds = 5),
            measure = msr("surv.cindex"),
            method_at = tnr("grid_search", resolution = 20, batch_size = 10),
            method_afs = NULL,
            fselector = FALSE,
            term_evals = NULL,
            workers = 10,
            outDir = NULL,
            parallel = TRUE,
            seed = 1993
        )


# Predictions
# =====
print("Making predictions in test data ...")
preds <- res$learner$predict_newdata(test)

normalize <- function(x, na.rm = TRUE) {
    return((x - min(x)) / (max(x) - min(x)))
}

scores <- data.frame(
    SampleID = rownames(test),
    Score = normalize(exp(preds$lp))
)

probs <- scores$Score
time <- test$Event_time
event <- test$Event

cindex <- Hmisc::rcorr.cens(probs, Surv(time, event), outx = FALSE)
print(paste0("C-Index in test set model is: ", 1 - cindex[1]))

write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))
