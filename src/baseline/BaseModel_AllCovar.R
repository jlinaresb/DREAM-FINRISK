# Copyright: Aki Havulinna @ THL, modified by Pande Putu Erawijantari
rm(list=ls())

library(survival)
#library(ResourceSelection)
#library(GGally)
#library(StepReg)
library(dplyr) # %>%
library(microbiome)
library(Hmisc) 
source("./modelperf_aki.R")

args=(commandArgs(TRUE))
PARAM <- list()
#this argument is to specify the path of input folder
#the input folder structure is similar to DreamHF.zip both for synthetic and real dataset
PARAM$folder.R <- paste0(args[1]) 
#one output file (score.csv) has to be created in Team_Name_Submission_Number folder
#in bellow example your submission name is : FINRISK_TEST_1, please change to your submission name
#please avoid using (.) in your team name
dir.create(file.path(PARAM$folder.R, "FINRISK_TEST_2","output"))
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
PARAM$folder.result <- paste0(PARAM$folder.data, "FINRISK_TEST_2/output/")

#functions
#dataframe filtering and modifications
transmerge <- function(O,S,core_taxa,fakenames=FALSE){ 
  O.new <- phyloseq(otu_table(as.matrix(O), taxa_are_rows=T))
  O.clr <- transform(O.new,"clr")
  O.clr <- prune_taxa(core_taxa, O.clr) # this is now not filtered to core taxa 
  if(fakenames==TRUE) {
    # change rownames with fake ids for convenience
    rownames(O.clr) <- paste0('bacteria',1:nrow(O.clr))
  }
  df <- cbind(meta(S), t(O.clr))
  df$Event_time<- as.numeric(df$Event_time)
  # exclude any minus or missing value
  df <- subset(df, Event_time > 0 & Event_time < 17)
  df <- subset(df, !is.na(Event_time))
  df <- subset(df, !is.na(Event)) 
  df = df %>% dplyr::select(!PrevalentHFAIL) 
  #replace NA with 0
  df <- df %>% replace(is.na(.), 0)
}


# load data
print("Load data")
S.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/pheno_test.csv"), 
                   row.names=1,)
S.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/pheno_training.csv"),
                   row.names=1,)

O.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/readcounts_test.csv"), 
                   row.names=1,)
O.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/readcounts_training.csv"), 
                    row.names=1,)


##define_taxa to filter
O.new <- phyloseq(otu_table(as.matrix(O.train), taxa_are_rows=T))
core_taxa <- core(O.new%>% transform("compositional"),detection = .1/100, prevalence = 1/10)%>% taxa_names()

df.train <- transmerge(O.train,S.train,core_taxa=core_taxa,fakenames=TRUE)
df.test <- transmerge(O.test,S.test,core_taxa=core_taxa,fakenames=TRUE)
predictors <- names(df.train)
predictors <- setdiff(predictors,c("PrevalentHFAIL","Event_time", "Event"))
basepredictors <- intersect(predictors,names(S.train))

survo <- 'Surv(Event_time,Event)~'
basemodel <- as.formula(paste0(survo,paste(basepredictors,collapse='+')))
# predicted absolute risk with training data
# the model, has to be calibrated and have decent Cindex value, i.e. be predictive
basemodel.fit.train <- coxph(basemodel,data=df.train)

basemodel.pred.test.obj <- Coxar.pred(basemodel.fit.train,newdata=df.test,years=15)  
#score file to test the scoring flow
res <- data.frame(SampleID=rownames(basemodel.pred.test.obj$y), 
                  Score=(basemodel.pred.test.obj$pred))  

write.csv(res, file=paste0(PARAM$folder.result, 
                           "scores.csv"), 
          quote=FALSE, row.names=FALSE)
