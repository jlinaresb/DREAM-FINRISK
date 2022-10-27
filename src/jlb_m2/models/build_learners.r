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