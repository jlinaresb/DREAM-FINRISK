# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_Final_Submission"
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
# source(file.path(code_path, "utils/predRes_helper.r"))
# source(file.path(code_path, "utils/holmes_test.r"))
# source(file.path(code_path, "utils/modelperf_aki.R"))


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

# Build cox with all train data
# ====
selected_betas <- models$lasso
train <- train[, match(
  c("Event", "Event_time", names(selected_betas)),
  colnames(train))]

surv <- "Surv(Event_time, Event) ~"
f <- as.formula(paste0(surv, paste(names(selected_betas),collapse='+')))
model <- coxph(f, train)

# Risk Prediction
# ====
risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = time)$surv)
}

pred <- risk(model, test, time = 15)

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)


# Save scores
scores <- data.frame(
    SampleID = rownames(test),
    Score = pred
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))