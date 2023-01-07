xgboost_surv_pipeline <- function(data,
                                 dataname,
                                 time,
                                 event,
                                 positive,
                                 removeConstant,
                                 normalize,
                                 filterFeatures,
                                 inner,
                                 measure,
                                 method_at,
                                 method_afs,
                                 term_evals,
                                 fselector,
                                 workers,
                                 outDir,
                                 parallel,
                                 seed) {
  set.seed(seed)
  # Make task
  task <- making_surv_task(data,
                           dataname,
                           time,
                           event)
  # Preprocess
  task <- preprocess(task,
                     removeConstant,
                     normalize,
                     filterFeatures)
  # Learner
  learner <- xgboost_surv(inner,
                         measure,
                         method_at,
                         method_afs,
                         term_evals,
                         fselector)
  # Parallelization
  if (parallel == TRUE) {
          future::plan(list(
          future::tweak("multisession", workers = 10))) # inner
  }
  # Autotuner
  rr <- learner$train(task)

  return(rr)
}