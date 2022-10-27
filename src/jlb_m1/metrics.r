setwd(here::here())
source("src/utils/hosmerTest.R")
setwd("src/jlb_m1/results/")

experiment <- readRDS("model_diff_exp.rds")

i <- 8
# Print C-Index
t_train <- experiment$train$Event_time
e_train <- experiment$train$Event
s_train <- experiment$prediction$scores_train[, i]

t_test <- experiment$test$Event_time
e_test <- experiment$test$Event
s_test <- experiment$prediction$scores_extval[, i]

# Concordance
require(survival)
require(Hmisc)
c_train <- rcorr.cens(exp(-s_train), Surv(t_train, e_train), outx = FALSE)
c_test <- rcorr.cens(exp(-s_test), Surv(t_test, e_test), outx = FALSE)

paste0("C-Index in train set model ",
        names(experiment$models)[i],
        " is: ", c_train[1])
paste0("C-Index in test set model ",
        names(experiment$models)[i],
        " is: ", c_test[1])
