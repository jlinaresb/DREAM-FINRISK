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

source("/utils/requirements.r")
source("/utils/importPseq.r")
source("/utils/prepro_functions.r")
source("/utils/co-abundances.r")
source("/utils/get_scores.r")
source("/utils/fit_model.r")
source("/utils/predRes_helper.r")

# Preprocess
# ======
source("/utils/preprocessing.r")

# Agglomerate by Species
train <- tax_glom(train, taxrank =  "Species")
train <- subset_taxa(train, Species != "s__")
test <- tax_glom(test, taxrank = "Species")
test <- subset_taxa(test, Species != "s__")

# Check data
# ======
stopifnot(ncol(x_train) == ncol(x_test))
stopifnot(colnames(x_train) == colnames(x_test))
stopifnot(nrow(x_train) == nrow(y_train))
stopifnot(nrow(x_test) == nrow(y_test))

# Creating train and test data
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

models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = "lasso")


environment(predRes2) <- asNamespace("biospear")
assignInNamespace("predRes", predRes2, ns = "biospear")
prediction <- predRes2(res = models,
                    method = "lasso",
                    traindata = train,
                    newdata = test,
                    int.cv = FALSE,
                    int.cv.nfold = 5,
                    time = seq(1, 16, 1),
                    trace = TRUE,
                    ncores = 20)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

risk <- exp(prediction$scores_extval)[, 1]
risk <- risk / risk + 1

# Save scores
scores <- data.frame(
    SampleID = rownames(test),
    Score = risk
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))