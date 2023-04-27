# Argument parsing
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2"
inputdir <- args[1]
#inputdir <- "/home/joselinares/git/DREAM-FINRISK/input_real/"
outputdir <- file.path(
  dirname(inputdir),
  "FINRISKDREAM_post1",
  "TEAM",
  Team_Name_Submission_Number,
  "output_refining_v1")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)

start <- Sys.time()

code_path <- "/mnt"
#code_path <- "/home/joselinares/git/DREAM-FINRISK/src/refining"

source(file.path(code_path, "utils/requirements.r"))
source(file.path(code_path, "utils/importPseq.r"))
source(file.path(code_path, "utils/prepro_functions.r"))
source(file.path(code_path, "utils/co-abundances.r"))
source(file.path(code_path, "utils/get_scores.r"))
source(file.path(code_path, "utils/fit_model.r"))

print("Start preprocessing...")

# Preprocess
# ======
source(file.path(code_path, "utils/preprocessing.r"))



# Fit model
# ========
method = c("lasso", "alassoL", "alassoU", "enet")
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

saveRDS(models, file = file.path(outputdir, "models.rds"))

# Build cox with all train data
# ====
risk <- function(model, newdata, time) {
  as.numeric(1-summary(
    survfit(model, newdata = newdata, se.fit = F, conf.int = F), 
    times = time, extend = TRUE)$surv)
}

for(i in seq_along(models)) {

  selected_betas <- models[[i]]
  model_name <- names(models)[i]

  train_ <- train[, match(
    c("Event", "Event_time", names(selected_betas)),
    colnames(train)
  )]

  surv <- "Surv(Event_time, Event) ~"
  f <- as.formula(paste0(surv, paste(names(selected_betas),collapse='+')))
  model <- coxph(f, train_)

  # Predict risk
  # ===
  pred <- risk(model, test, time = 15)

  # Save scores
  scores <- data.frame(
      SampleID = rownames(test),
      Score = pred
  )
  write.csv(scores, quote = FALSE, row.names = FALSE,
            file = file.path(
              outputdir, 
              paste0("scores_", model_name, ".csv")))

  print(paste0("Predicting with ", model_name))
}

end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)
