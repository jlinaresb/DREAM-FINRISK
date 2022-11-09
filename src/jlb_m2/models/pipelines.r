setwd(here::here())
source("src/jlb_m2/models/build_learners.r")
source("src/jlb_m2/models/pipeline_utils.r")

rf_pipeline <- function(data,
                        dataname,
                        target,
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
  task <- making_task(data,
                      dataname,
                      target,
                      positive)
  # Preprocess
  task <- preprocess(task,
                     removeConstant,
                     normalize,
                     filterFeatures)
  # Learner
  learner <- randomForest(inner,
                          measure,
                          method_at,
                          method_afs,
                          term_evals,
                          fselector)
  # Parallelization
  if (parallel == TRUE) {
        future::plan(list(
            future::tweak("multisession", workers = 20))) # inner
  }
  # Autotuner
  model <- learner$train(task)
  saveRDS(model,
          file = paste0(outDir, "/model_randomForest_", dataname, ".rds"))

  return(model)
}


surv_pipeline <- function(data,
                        dataname,
                        time,
                        event,
                        positive,
                        removeConstant,
                        normalize,
                        filterFeatures,
                        type,
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
  learner <- ann_surv(type,
                      inner,
                      measure,
                      method_at,
                      method_afs,
                      term_evals,
                      fselector)
  # Parallelization
  if (parallel == TRUE) {
        future::plan(list(
            future::tweak("multisession", workers = 20))) # inner
  }
  # Autotuner
  model <- learner$train(task)
  saveRDS(model,
          file = paste0(outDir, "/model_randomForest_", dataname, ".rds"))

  return(model)
}
