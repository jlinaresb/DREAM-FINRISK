# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_4"
inputdir <- args[1]
outputdir <- file.path(inputdir, Team_Name_Submission_Number, "output")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)

start <- Sys.time()

code_path <- "/mnt"

source(file.path(code_path, "utils/requirements.r"))
source(file.path(code_path, "utils/importPseq.r"))
source(file.path(code_path, "utils/prepro_functions.r"))
source(file.path(code_path, "utils/co-abundances.r"))
source(file.path(code_path, "utils/get_scores.r"))
source(file.path(code_path, "utils/fit_model.r"))
source(file.path(code_path, "utils/predRes_helper.r"))
source(file.path(code_path, "utils/holmes_test.r"))

# Preprocess
# ======
source("/mnt/utils/preprocessing.r")


# Fit model
# ========
method = "lasso"
cvrts <- c("Age", "BodyMassIndex",
           "SystolicBP", "NonHDLcholesterol",
           "Sex", "disbiosis")

print("Fitting the model ...")
models <- fit_biospear(data = train,
                       biomarkers = setdiff(colnames(train),
                                            c(cvrts, colnames(y_train))),
                       surv = c("Event_time", "Event"),
                       cvrts = cvrts,
                       inter = FALSE,
                       methods = method)

print(models)

environment(predRes2) <- asNamespace("biospear")
assignInNamespace("predRes", predRes2, ns = "biospear")
print("Predict in test data ...")
prediction <- predRes2(res = models,
                    method = method,
                    traindata = train,
                    newdata = test,
                    int.cv = FALSE,
                    int.cv.nfold = 5,
                    time = seq(1, 15, 1),
                    trace = TRUE,
                    ncores = 8)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

normalize <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}

pred <- as.vector(prediction)
pred <- normalize(pred)

# C-Index in test
print(rcorr.cens(-pred, Surv(test$Event_time, test$Event))[1])
# Hoslem test
print(HosLem.test(Surv(test$Event_time, test$Event), pred)$pval)


# Save scores
scores <- data.frame(
    SampleID = rownames(test),
    Score = pred
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))