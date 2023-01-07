# Argument parsing
start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_3"
inputdir <- args[1]
outputdir <- file.path(inputdir, Team_Name_Submission_Number, "output")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)
code_path <- "/mnt"

source(file.path(code_path , "utils/requirements.r"))
source(file.path(code_path , "utils/importPseq.r"))
source(file.path(code_path , "utils/get_scores.r"))
source(file.path(code_path , "utils/co-abundances.r"))
source(file.path(code_path , "utils/holmes_test.r"))

# Preprocess
# =======
source(file.path(code_path , "utils/preprocessing.r"))
#load("~/git/DREAM-FINRISK/tmp/data_new.RData")

# Run Genetic Algorithm
# ========
# Split train
set.seed(1993)
train.index <- createDataPartition(train$Event, p = .7, list = FALSE)
trainPart <- train[train.index, ]
testPart  <- train[-train.index, ]

nbits <- ncol(train) - 2
fitness <- function(col) {
  
  # Create surv data
  surv_train <- Surv(time = trainPart[, "Event_time"],
                     event = trainPart[, "Event"])
  surv_test  <- Surv(time = testPart[, "Event_time"],
                     event = testPart[, "Event"])
  
  # Subset columns
  train_filt <- subset(trainPart, select = -c(Event_time, Event))[as.logical(col)]
  
  # get columns desired
  variables <- colnames(train_filt)
  
  # Formula
  f <- as.formula(paste0("surv_train ~ ", paste0(variables, collapse = " + ")))

  # Fit model
  fit <- coxph(f, trainPart)
  
  # Prediction in partition test
  prediction <- predict(fit, testPart, type = "risk")
  
  # Scale predictions
  prediction <- sigmoid(prediction, method = "logistic")
  
  # Harrel-C
  c <- Hmisc::rcorr.cens(-prediction, surv_test)[1]
  
  # Hoslem test
  h <- HosLem.test(surv_test, prediction)$pval
  
  # Add penalty in case of h <<< 0.05 (modify?)
  if (h < 1e-150) {
    c <- c - 0.4
  } else if (h < 1e-100) {
    c <- c - 0.3
  } else if (h < 1e-50) {
    c <- c - 0.2
  } else if (h > 1e-50 & h < 1e-10) {
    c <- c - 0.1
  } else if (h > 1e-10 & h < 0.05) {
    c <- c + 0.05
  } else  if (h > 0.05) {
    c <- c + 0.1
  }
  
  return(h)
}

# Genetic algorithm
print("Run Genetic Algorithm ...")
result <- ga(
  type = "binary",
  fitness = fitness,
  nBits = nbits,
  popSize = 200,
  pmutation = 0.5,
  pcrossover = 0.9,
  run = 10,
  elitism = 10, 
  maxiter = 1000,
  monitor = TRUE,
  parallel = parallel::detectCores(),
  seed = 1993)

saveRDS(result, file = "~/git/DREAM-FINRISK/tmp/ga_cox_m3.rds")

# Extract best solution
solution <- which.min(rowSums(result@solution))
solution <- result@solution[solution, ]

# Retrain with all train subset
surv_train <- Surv(time = train[,"Event_time"], event = train[,"Event"])
surv_test <- Surv(time = test[,"Event_time"], event = test[,"Event"])

x_train <- subset(train, select = -c(Event_time, Event))
x_test <- subset(test, select = -c(Event_time, Event))

x_train <- x_train[, as.logical(solution), drop=F]

f <- as.formula(paste0("surv_train ~ ", paste0(names(x_train), collapse=" + ")))

# Fit the model
fit <- coxph(f, train)

# Prediction in partition test
prediction <- predict(fit, test, type = "lp")

# Scale predictions
prediction <- sigmoid(prediction, method = "logistic")

# C-Index in test
c <- rcorr.cens(-prediction, surv_test)
print(paste0("C-Index in test set: ", c[1]))

# Holmes test
h <- HosLem.test(surv_test, prediction)
print(paste0("Holmes p-value in test set: ", h$pval))

# Ellapsed time
end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

# Save scores
scores <- data.frame(
  SampleID = rownames(test),
  Score = prediction,
  row.names = 1:nrow(test)
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))
