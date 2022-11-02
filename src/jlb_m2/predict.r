setwd(here::here())
source("src/jlb_m2/requirements.r")
model <- readRDS("src/jlb_m2/results/model_randomForest_diff_rf.rds")
rr <- readRDS("src/jlb_m2/results/rsmp_randomForest_diff_rf.rds")


# In resample results
model_x <- rr$result$learners[[1]]$learner
model_x$predict_newdata(data_test, task = NULL)

# In autotuner
pred <- model$learner$predict_newdata(data_test, task = NULL)


names(model)
source("src/utils/importPseq.r")
source("src/jlb_m2/models/pipeline_utils.r")
test <- pseq(subset = "test")
names(rr$learners[[1]])
class(rr$learner)
task <- TaskClassif$new(id = "test",
                    backend = data_test,
                    target = "Event",
                    positive = NULL)




rr$learners[[1]]$train(model$task)
