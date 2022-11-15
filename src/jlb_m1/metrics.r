require(survival)
require(Hmisc)
setwd(here::here())
source("src/utils/hosmerTest.R")
setwd("src/jlb_m1/results/")

experiment <- readRDS("model_pca.rds")
t_train <- experiment$train$Event_time
e_train <- experiment$train$Event

t_test <- experiment$test$Event_time
e_test <- experiment$test$Event

m <- names(experiment$models)
lapply(m, function(i) {
    s_train <- experiment$prediction$scores_train[, i]
    s_test <- experiment$prediction$scores_extval[, i]
    # Concordance
    c_train <- rcorr.cens(exp(-s_train), Surv(t_train, e_train), outx = FALSE)
    c_test <- rcorr.cens(exp(-s_test), Surv(t_test, e_test), outx = FALSE)
    # Print C-Index
    print(paste0("C-Index in train set model ", i, " is: ", c_train[1]))
    print(paste0("C-Index in test set model ", i, " is: ", c_test[1]))
})

experiment$models
