setwd(here:here())
source("src/jlb_m2/models/build_learners.r")
source("src/jlb_m2/models/pipelines_utils.r")

rf_pipeline <- function(data,
                        dataname,
                        target,
                        positive,
                        removeConstant,
                        normalize,
                        filterFeatures,
                        inner,
                        outer,
                        measure,
                        method_at,
                        method_afs,
                        term_evals,
                        fselector,
                        workers,
                        outDir,
                        parallel,
                        folds,
                        batch_size,
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
            future::tweak("multisession", workers = folds),  # outer
            future::tweak("multisession", workers = batch_size))) # inner
  }
  # Resampling
  rr <- resample(task,
                 learner,
                 resampling = outer,
                 store_models = TRUE)
  # Save resampling object
  res <- list(task = task,
              result = rr)
  saveRDS(res,
          file = paste0(outDir, "/rsmp_randomForest_", dataname, ".rds"))
}