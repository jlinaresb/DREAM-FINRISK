setwd("~/git/DREAM-FINRISK/src/jlb_m1/results/")
experiment <- readRDS("model_diff_exp.rds")

i <- 1
# Print C-Index
t_train <- experiment$train$Event_time
e_train <- experiment$train$Event
s_train <- experiment$prediction$scores_train[, i]

t_test <- experiment$test$Event_time
e_test <- experiment$test$Event
s_test <- experiment$prediction$scores_extval[, i]

cindex_train <- survcomp::concordance.index(s_train, t_train, e_train)
cindex_train$c.index
cindex_test <- survcomp::concordance.index(s_test, t_test, e_test)
cindex_test$c.index


# Concordance
require(survival)
require(Hmisc)
xx <- rcorr.cens(exp(-s_test), Surv(t_test, e_test), outx = FALSE)
xx[1]
