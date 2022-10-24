# load libraries 
library(tidyverse)
library(survminer)
library(survival)
library(microbiome)
library(magrittr)
library(randomForestSRC)
set.seed(198657)

mainDir <- getwd()
subDir <- "output"

# create output folder
if (file.exists(subDir)){
  setwd(file.path(mainDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir))
}
# set parameters
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "output/")
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
source(paste0(mainDir,"/src/hosmerTest.R"))

# load data
print("Load data")

S.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/pheno_training.csv"),
                   row.names=1,)
S.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/pheno_test.csv"), 
                   row.names=1,)
O.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/readcounts_training.csv"), 
                    row.names=1,)
O.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/readcounts_test.csv"), 
                   row.names=1,)
###############################################################################
# set parameters
alpha_level = 1
endpoints <- c("Event_time", "Event")

# change rownames with fake ids for convenience
rows <- rownames(O.train)  
seq=seq(1:c(length(rows)))
fakename<- sapply(seq,
                  function(x) paste('bacteria', x))
fakename = gsub(" ","",as.character(fakename))
rownames(O.train) <- fakename

df.train <- cbind(meta(S.train), t(O.train))
df.train$Event_time<- as.numeric(df.train$Event_time)
# exclude any minus value
df.train <- subset(df.train, Event_time > 0 & Event_time < 17)
df.train <- subset(df.train, !is.na(Event_time))
df.train <- subset(df.train, !is.na(Event))  

# define potential predictors
predictors <- c(colnames(S.train), rownames(O.train))

# choose only HF cases and eliminate "NA"
df.train <- subset(df.train, select=predictors, PrevalentHFAIL==0&!is.na(Event))

# remove PrevalentHFAIL from all matrix
predictors <- predictors[! predictors %in% "PrevalentHFAIL"]
df.train = df.train %>% 
  select(!PrevalentHFAIL) 
###############################################################################
# Fix and combine train and test datasets
# change rownames with fake ids for convenience
rownames(O.test) <- fakename

df.test <- cbind(meta(S.test), t(O.test))
df.test$Event_time <- as.numeric(df.test$Event_time)
# exclude any minus value
df.test <- subset(df.test, Event_time > 0 & Event_time < 17)
df.test = df.test %>% 
  select(!PrevalentHFAIL) 

# remove unnecessary files and save the main files
rm(O.train, O.test, S.test, S.train,fakename )
#save(df.train, df.test, file=paste0(PARAM$folder.result,"dream_synthetic_files.RData"))

###############################################################################
# Univariate Cox regression analysis*******************************************

# A Cox regression of time to death on the time-constant covariates is specified as follow:
# cox ph model: allows both categorical and numeric variables 
# with categorical variable
# we cant estimate survival/hazard  but estimate hazard ratio
# the baseline survival function

predictors <- setdiff(predictors, endpoints)

cox.mod.sex <- coxph(Surv(Event_time, Event) ~ Sex+Age,
                     data=df.train) 

summary(cox.mod.sex)
print(cox.zph(cox.mod.sex))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.subset.pdf"))
ggforest(cox.mod.sex, data=df.train)
dev.off() # Turn the PDF device off

# All meta variables *****************************
cox.mod.total.meta <- coxph(Surv(Event_time, Event) ~ BodyMassIndex + Smoking + 
                              PrevalentDiabetes + SystolicBP + BPTreatment + 
                              NonHDLcholesterol + PrevalentCHD + Sex  + Age, 
                            data=df.train)
summary(cox.mod.total.meta)
print(cox.zph(cox.mod.total.meta))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.total.pdf"))
ggforest(cox.mod.total.meta, data=df.train)
dev.off() 

save(cox.mod.sex, cox.mod.total.meta, 
     file=paste0(PARAM$folder.result, "cox_models.RData"))

# Multivariate models **********************************************************
# Calculate HarrelC
model.list <- list(cox_subset = cox.mod.sex, 
                   cox_meta_total = cox.mod.total.meta)

collect.harrelc <- data.frame()
collect.hoslem <- data.frame()

for (model in model.list){  
  print(model$formula)
  scores=as.data.frame(predict(model, 
                               newdata=df.test,     
                               se.type="expected")) 
  scores = scores %>%
    set_colnames(c("Score")) %>%
    rownames_to_column(var = "SampleID")
  labels <- df.test
  labels$SampleID<- rownames(df.test)
  
  # Align the user provided scores with the true event times
  true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
  
  # Calculate Harrell's C statistics
  HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
  collect.harrelc <- rbind(collect.harrelc, HarrellC["C Index"])
  colnames(collect.harrelc) <- "HarrellC"

  # Calculate Hoslem
  pred <- Coxar(model, years=10, bh=NULL)
  test <- HosLem.test(model$y, pred, plot=FALSE)
  hoslem_p <-  test$pval
  collect.hoslem <- rbind(collect.hoslem, hoslem_p)
  colnames(collect.hoslem) <- "hoslem_p"
}

rownames(collect.hoslem) <- names(model.list)
collect.hoslem$hoslem_p.adj <- p.adjust(collect.hoslem$hoslem_p)
rownames(collect.harrelc) <- names(model.list)
totalscores <-cbind(collect.harrelc, collect.hoslem)

# export file
write.csv(totalscores, file=paste0(PARAM$folder.result, 
                                    Sys.Date(),"_",
                                    "cox_scores.csv"))

# plot cox 
# Create the new df for sex  
sex_df <- with(df.train,
               data.frame(Sex = c(0, 1), 
                          Age = rep(mean(Age, na.rm = TRUE), 2)))

fit <- survfit(cox.mod.sex, newdata = sex_df)
plot.agesex=ggsurvplot(fit, conf.int = FALSE, data = sex_df, palette = "Dark2", 
                 censor = FALSE, surv.median.line = "hv")

ggsave(file = paste0(PARAM$folder.result, 
                     Sys.Date(), "_",
                     "ggsurvplot_byAge.pdf"), plot.agesex$plot)

###############################################################################
# Random survival forest********************************************************
set.seed(9876)

# Labels are given for parameter optimization
dff1 <- as.data.frame(df.train[, c(endpoints, predictors)])
rownames(dff1)<-NULL

# Labels are hidden for the second data subset 
dff2 <- as.data.frame(df.test[, predictors])
rownames(dff2)<-NULL

# Train Random Survival Model with the labeled training data
rf1 <- rfsrc(Surv(Event_time, Event) ~ .,
             data = dff1,
             importance = FALSE)

# Predictions for the test data
samples <- which(rowSums(is.na(dff2))==0)
pred2 <- predict(rf1, dff2[samples, ])
scores <- pred2$predicted; names(scores) <- rownames(df.test)[samples]
res <- data.frame(SampleID=rownames(df.test), 
                  Score=sample(scores[rownames(df.test)]))  # Fully random
res <- res %>%  drop_na()
# Write the scoring table for evaluation
write.csv(res, file=paste0(PARAM$folder.result, 
                           Sys.Date(), 
                           "random_forest_scores.csv"), 
          quote=FALSE, row.names=FALSE)

# Harrells C *******************************************************************
labels=df.train
labels$SampleID<- rownames(df.train)

# Align the user provided scores with the true event times
true.scores <- as.numeric(labels[res$SampleID,"Event_time"])

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(res$Score, true.scores, outx=FALSE)
print(C)

Cindex.train <- get.cindex(rf1$yvar[,1], rf1$yvar[,2], rf1$predicted.oob)

# Hoslem test ***************************************************************** 
pred <- Coxar(model, years=10, bh=NULL)
test <- HosLem.test(model$y, pred, plot=FALSE)
# subsample
oo <- subsample(rf1)

# combine all
stats_rfs <- data.frame(harrellc.test = C[[1]], 
                        harrellc.train = Cindex.train,
                        hosLem.pval = test$pval) 
# export stats
write.csv(stats_rfs, file=paste0(PARAM$folder.result, 
                         Sys.Date(), 
                         "random_forest_stats.csv"))

# plot
pdf(paste0(PARAM$folder.result, "RSF_overview.pdf"), 
    width = 10, 
    height = 8)
par(oma = c(0.5, 0.5, 0.5, 0.5))
plot.survival(rf1, cens.model = "rfsrc",  show.plots = TRUE)
dev.off()

pdf(paste0(PARAM$folder.result, "Variable Importance RSF.pdf"), 
    width = 10, 
    height = 15)
par(oma = c(0.5, 0.5, 0.5, 0.5))
par(cex.axis = 2.0,
    cex.lab = 2.0,
    cex.main = 2.0,
    mar = c(6.0,17,1,1),
    mgp = c(4, 1, 0))
plot.subsample(oo, alpha = .01, standardize = TRUE, normal = TRUE, pmax = 20)
dev.off()


