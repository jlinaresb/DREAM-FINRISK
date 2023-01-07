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
    max_depth = p_int(lower = 3, upper = 16)
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