require(survival)
require(Hmisc)
setwd(here::here())
source("src/utils/hosmerTest.R")
setwd("src/jlb_m1/results/")

load("src/jlb_m1/results/res_coab_lasso.rds")

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
