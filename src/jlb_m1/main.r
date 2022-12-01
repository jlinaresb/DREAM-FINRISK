# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2"
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
source("/mnt/utils/co-abundances.r")
source("/mnt/utils/get_scores.r")
source("/mnt/utils/fit_model.r")
source("/mnt/utils/predRes_helper.r")

# Preprocess
# ======
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

# Fit model
# ========
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex")

print("Fitting the model ...")
models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = "lasso")

print(models)

environment(predRes2) <- asNamespace("biospear")
assignInNamespace("predRes", predRes2, ns = "biospear")
print("Predict in test data ...")
prediction <- predRes2(res = models,
                    method = "lasso",
                    traindata = train,
                    newdata = test,
                    int.cv = FALSE,
                    int.cv.nfold = 5,
                    time = seq(1, 15, 1),
                    trace = TRUE,
                    ncores = 20)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

risk <- exp(prediction$scores_extval)[, 1]
normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}
risk <- normalize(risk)

# Save scores
scores <- data.frame(
    SampleID = rownames(test),
    Score = risk
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))