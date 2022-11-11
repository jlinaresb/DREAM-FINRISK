source("src/jlb_m2/models/tuning_utils.r")

# Random Forest
# ====
randomForest <- function(inner,
                         measure,
                         method_at,
                         method_afs,
                         term_evals,
                         fselector) {
  # Make learner
  learner <- lrn("classif.randomForest",
                  predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    mtry = p_int(lower = 3, upper = 15),
    nodesize = p_int(lower = 1, upper = 5),
    ntree = p_int(lower = 990, upper = 1000)
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}



ann_surv <- function(type,
                     inner,
                     measure,
                     method_at,
                     method_afs,
                     term_evals,
                     fselector) {

  # Make learner
  learner <- lrn(type,
                 epochs = 100,
                 optimizer = "adam")

  # Hyperparameter space
  ps <- ps(
    dropout = p_dbl(lower = 0, upper = 1),
    weight_decay = p_dbl(lower = 0, upper = 0.5),
    learning_rate = p_dbl(lower = 0, upper = 1),
    nodes = p_int(lower = 1, upper = 32),
    k = p_int(lower = 1, upper = 20),
    batch_size = p_int(lower = 24, upper = 500)
  )
  ps$trafo <- function(x, param_set) {
    x$num_nodes <- rep(x$nodes, x$k)
    x$nodes <- x$k <- NULL
    return(x)
  }
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    method_at,
                    fselector,
                    method_afs)

  return(afs)
}


xgboost_surv <- function(inner,
                         measure,
                         method_at,
                         method_afs,
                         term_evals,
                         fselector) {

  # Make learner
  learner <- lrn("surv.xgboost",
                 nthread = 10)

  # Hyperparameter space
  ps <- ps(
    eta = p_dbl(lower = 0.01, upper = 0.3),
    gamma = p_dbl(lower = 0, upper = 9),
    colsample_bytree = p_dbl(lower = 0.5, upper = 1),
    max_depth = p_int(lower = 3, upper = 10)
  )

  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    method_at,
                    fselector,
                    method_afs)

  return(afs)
}