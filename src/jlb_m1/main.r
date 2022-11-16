setwd(here::here())
experiment_id <- "coab_lasso"
outdir <- "src/jlb_m1/results/"
save <- FALSE

start <- Sys.time()
source("src/jlb_m1/requirements.r")
source("src/utils/importPseq.r")
source("src/utils/prepro_functions.r")
source("src/utils/co-abundances.r")
source("src/utils/get_scores.r")
source("src/jlb_m1/fit_model.r")
source("src/jlb_m1/predRes_helper.r")

# Preprocess
# ======
source("src/preprocessing/preprocessing.r")

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
test$Even_time[which(test$Event_time < 0), ] <- 15

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

if (save == TRUE)  {
    save(train,  test, models, prediction,
         file = paste0("src/jlb_m1/results/res_", experiment_id, ".rds"))
}

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

# Save scores
scores <- data.frame(
    SampleID = rownames(test),
    Score = exp(-prediction$scores_extval)[, 1]
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = "output/scores.csv")