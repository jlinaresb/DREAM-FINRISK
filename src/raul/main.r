# Argument parsing
# args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_3"
# inputdir <- args[1]
inputdir <- "~/git/DREAM-FINRISK/"
outputdir <- file.path(inputdir, Team_Name_Submission_Number, "output")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)

start <- Sys.time()

code_path <- "~/git/DREAM-FINRISK/src/raul/"

source(file.path(code_path, "requirements.r"))
source(file.path(code_path, "importPseq.r"))
source(file.path(code_path, "prepro_functions.r"))
source(file.path(code_path, "co-abundances.r"))
source(file.path(code_path, "get_scores.r"))

# Preprocess
# ======
source("preprocessing.r")


# Check data
# ======
stopifnot(ncol(x_train) == ncol(x_test))
stopifnot(colnames(x_train) == colnames(x_test))
stopifnot(nrow(x_train) == nrow(y_train))
stopifnot(nrow(x_test) == nrow(y_test))

# Creating train and test data
print("Creating train and test data ...")
train <- cbind.data.frame(x_train, y_train)
test  <- cbind.data.frame(x_test, y_test)

# What do we do with ...
# ======
# Negatives survival values in train and test?
train <- train[-which(train$Event_time < 0), ]
test$Event_time[which(test$Event_time < 0)] <- 15

# NA's values in train and test?
train <- train[complete.cases(train), ]
test <- missRanger(test, pmm.k = 10, seed = 153)

library(survival)
library(dynpred)
library(caret)
library(ResourceSelection)
# filter nas
train=na.omit(train)
test=na.omit(test)
# remove heartfail
train=train[,-c(7)]
test=test[,-c(7)]
testF=test
# create index for train partition
train.index <- createDataPartition(train$Event, p = .7, list = FALSE)

# fitness function for genetic alg
fitcox <- function(col){
  # split train into train and test
  trainPart <- train[ train.index,]
  testPart  <- train[-train.index,]
  
  # create surv data
  surv_train=survival::Surv(time=trainPart[,"Event_time"],event=trainPart[,"Event"])
  surv_test=survival::Surv(time=testPart[,"Event_time"],event=testPart[,"Event"])
  
  # get columns desired
  variables=colnames(subset(trainPart,select=-c(Event_time,Event)))[as.logical(col)]
  
  form=paste0("surv_train ~ ", paste0(variables, collapse=" + "))
  #form=gsub('"','',form)
  #form=gsub(' ','',form)
  # desing formula
  formulacox <- formula(form)
  #formulacox <- formula("surv_train ~ .")
  # define cox model
  coxmodel <- coxph(formulacox,subset(trainPart,select=-c(Event_time,Event)))
  # predict test with cox model
  predicted <- predict(coxmodel,subset(testPart,select=-c(Event_time,Event)),type = "risk")
  
  # compute cindex
  #cindex <- Hmisc::rcorr.cens(-predicted,surv_test, outx = FALSE)[1]
  cindex=1-Hmisc::rcorr.cens(x = predicted,surv_test)[1]
  hoslem_test=hoslem.test(g = 10,x = predicted,y = testPart$Event)
  return(hoslem_test$p.value)
}

#start genetic alg
start_time <- Sys.time()
alg_get=(GA::ga(
  type = "binary", #tipo de datos
  fitness = fitcox, #funcion fitness
  nBits = 41, #numero de variables de entrada
  run = 50, #numero de iteraciones que el algoritmo espera a que mejore el resultado
  maxiter = 50, #numero de iteraciones
  #monitor = plot, #mostrar una grafica al final de cada iteracion
  parallel = 6, #permitir procesamiento en paralelo
  keepBest = T, #mantener la mejor solucion al final
  seed = 84221, #fijar la semilla para poder ser reproducible
  elitism = 15,
  names = colnames(subset(test,select=-c(Event_time,Event))),
  popSize = 1000, #numero de soluciones,
  pmutation = 0.2
))
end_time <- Sys.time()

# get best solution, last by default
bestsol=0
for(i in nrow(alg_get@solution)){
  if(ncol(alg_get@solution)>=sum(alg_get@solution[i,])){
    bestsol=alg_get@solution[i,]
  }
}
print("best solution")
print(names(bestsol)[as.logical(bestsol)])
# function to normalize prediction
normalize <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}

surv_train=survival::Surv(time=train[,"Event_time"],event=train[,"Event"])
surv_test=survival::Surv(time=testF[,"Event_time"],event=testF[,"Event"])

best_train=train[,as.logical(bestsol), drop=F]
variables=names(bestsol)[as.logical(bestsol)]

form=paste0("surv_train ~ ", paste0(variables, collapse=" + "))
#form=gsub('"','',form)
formulacox <- formula(form)
coxmodel <- coxph(formulacox,train[,variables])
predicted <- predict(coxmodel,subset(testF,select=-c(Event_time,Event)),type = "risk")
predictedNormalized <- normalize(predicted)
print(hoslem.test(g = 10,x = predictedNormalized,y = testF$Event))
# Save scores
scores <- data.frame(
  SampleID = rownames(test),
  Score = predictedNormalized
)
print(head(scores))
write.csv(scores, quote = FALSE, row.names = FALSE,
          file = file.path(outputdir, "scores.csv"))