# TODO:
# - Return scores.csv

setwd(here::here())
experiment_id <- "prueba"
outdir <- "src/jlb_m1/results/"

start <- Sys.time()
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

train <- cbind.data.frame(x_train, y_train)
test  <- cbind.data.frame(x_test, y_test)

# What do we do with the negatives?
# ======
# Removing them!
train <- train[-which(train$Event_time < 0), ]
test <- test[-which(test$Event_time < 0), ]

# NA"s??
train <- train[complete.cases(train), ]
test <- test[complete.cases(test), ]

# Removing PrevalentHFAIL (only 0)
train <- train[, -grep("PrevalentHFAIL", colnames(train))]
test <- test[, -grep("PrevalentHFAIL", colnames(test))]

save(train, test, file = "tmp/data_clusters.RData")

# Fit model
# ========
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex")
tt <- "BPTreatment"

methods <- c("alassoL", "alassoR",
            "alassoU", "enet", "gboost",
            "lasso", "lasso-1se",
            "lasso-AIC", "lasso-BIC", "lasso-HQIC",
            "lasso-pct", "lasso-pcvl", "lasso-RIC", "modCov",
            "PCAlasso", "PLSlasso", "ridge", "ridgelasso",
            "stabSel", "uniFDR")
models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(x_train), c(cvrts, tt)),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = TRUE,
                       treatment = tt,
                       methods = methods)

environment(predRes2) <- asNamespace("biospear")
assignInNamespace("predRes", predRes2, ns = "biospear")
prediction <- predRes2(res = models,
                    method = methods,
                    traindata = train,
                    newdata = test,
                    int.cv = FALSE,
                    int.cv.nfold = 5,
                    time = seq(1, 16, 1),
                    trace = TRUE,
                    ncores = 20)


t_train <- train$Event_time
e_train <- train$Event

t_test <- test$Event_time
e_test <- test$Event

require(Hmisc)
m <- names(models)
lapply(m, function(i) {
    s_train <- prediction$scores_train[, i]
    s_test <- prediction$scores_extval[, i]
    # Concordance
    c_train <- rcorr.cens(exp(-s_train), Surv(t_train, e_train), outx = FALSE)
    c_test <- rcorr.cens(exp(-s_test), Surv(t_test, e_test), outx = FALSE)
    # Print C-Index
    print(paste0("C-Index in train set model ", i, " is: ", c_train[1]))
    print(paste0("C-Index in test set model ", i, " is: ", c_test[1]))
})



end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

# Save data and models
# =======
saveRDS(list(
            train = train,
            test = test,
            models = models,
            prediction = prediction),
            file = paste0(outdir, "model_", experiment_id, ".rds"))

scores <- data.frame(
    SampleID = rownames(test),
    Score = exp(-prediction$precitions)
)

write.csv(scores, quote = FALSE, row.names = FALSE,
          file = "output/scores.csv")