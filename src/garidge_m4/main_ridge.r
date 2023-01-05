# Argument parsing
start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
Team_Name_Submission_Number <- "SB2_3"
# inputdir <- args[1]                UNCOMMENT!!!
inputdir <- "~/git/DREAM-FINRISK"  # DELETE!!
outputdir <- file.path(inputdir, Team_Name_Submission_Number, "output")

print(Team_Name_Submission_Number)
print(inputdir)
print(outputdir)

dir.create(outputdir, recursive = TRUE)
code_path <- "~/git/DREAM-FINRISK/src/rmoo/" # change to relative singularity path 

source(file.path(code_path , "utils/requirements.r"))
source(file.path(code_path , "utils/importPseq.r"))
source(file.path(code_path , "utils/get_scores.r"))
source(file.path(code_path , "utils/co-abundances.r"))

# Preprocess
# =======
source(file.path(code_path , "utils/preprocessing.r"))


# Clean ls
rm(list = setdiff(ls(), c("train", "test")))
gc()


# Run Genetic Algorithm
# ========
# Split train
set.seed(1993)
train.index <- createDataPartition(train$Event, p = .7, list = FALSE)

nbits <- ncol(train) - 2
fitness <- function(col) {
  
  trainPart <- train[train.index, ]
  testPart <- train[-train.index, ]
  
  # create surv data
  surv_train <- Surv(time = trainPart[, "Event_time"],
                     event = trainPart[, "Event"])
  surv_test <- Surv(time = testPart[, "Event_time"],
                    event = testPart[, "Event"])
  
  # get columns desired
  variables <- colnames(subset(trainPart, select = -c(Event_time, Event)))[as.logical(col)]
  
  # Formula
  # f <- as.formula(paste0("surv_train ~ ", paste0(variables, collapse = " + ")))

  # Fit model
  cvfit <- cv.glmnet(
    x = as.matrix(subset(trainPart, select=-c(Event_time, Event))),
    y = surv_train,
    type.measure = "C",
    family = "cox",
    alpha = 0)
  
  # Prediction in partition test
  prediction <- predict(cvfit,
                        as.matrix(subset(testPart, select=-c(Event_time, Event))),
                        s = "lambda.min")

  # Scale predictions
  normalize <- function(x, na.rm = TRUE) {
    return((x - min(x)) / (max(x) - min(x)))
  }
  prediction <- normalize(prediction)
  prediction <- as.vector(prediction)
  
  # Harrel-C
  c <- Hmisc::rcorr.cens(-prediction, surv_test)[1]
  
  # Hoslem test
  # h <- hoslem.test(x = testPart$Event, y = -prediction, g = 10)$p.value
  # Holmes test
  grpkm <- function(grps,surv,weights=NULL) {
    if(!is.Surv(surv)) stop('Please use a survival object')
    res <- NULL
    for(i in sort(unique(grps))) {
      grp <- grps==i
      kpla <- survfit(surv[grp,]~1,weights=weights[grp])
      last <- length(kpla$time)
      n <- sum(grp)
      events <- sum(kpla$n.event)
      KM <- kpla$surv[last]
      std.err <- kpla$stderr
      res <- rbind(res,cbind(n,events,KM,std.err))
    }
    res
  }
  tiles <- c('','','tertiles','quartiles','quintiles','sextiles','septiles','octiles','noniles','deciles')
  HosLem.test <- function(surv,pred,plot=TRUE,DF.reduce=2) {
    # Cook-Ridker version test uses DF.reduce=2 (default)
    # Cook NR, Ridker PM. (2009) Ann Intern Med. 150(11): 795-802
    # D'Agostino-Nam version uses DF.reduce=1
    # About the differences: http://hdl.handle.net/1773/22648
    # (Guffey D., Hosmer-Lemeshow goodness-of-fit test: Translations to the Cox Proportional Hazards Model)
    # if plot is a name, a jpg plot is saved with this name (.jpg will be added)
    if(!DF.reduce%in%c(1,2)) stop('Please specify DF.reduce = 1 or 2')
    if(!is.Surv(surv)) stop('Please use a survival object')
    version <- c('D\'Agostino-Nam','Cook-Ridker')[DF.reduce]
    grp <- as.numeric(cut2(pred,g=10))
    pj <- as.numeric(by(pred,grp,mean))
    idx <- TRUE
    while(any(idx)&length(pj)>3) {
      kmj <- 1-grpkm(grp,surv)[,3] # failure, not survival
      pj <- as.numeric(by(pred,grp,mean))
      nj <- table(grp)
      idx <- grp%in%which(nj<5|pj*nj<5|nj*kmj<5)&grp<max(grp)
      if(any(idx)) {
        grp[idx] <- grp[idx]+1
        grp <- match(grp,sort(unique(grp)))
      }
    }
    dan_chi2 <- sum(nj*(kmj-pj)^2/(pj*(1-pj)))
    pval <- pchisq(dan_chi2,length(pj)-DF.reduce,lower.tail=FALSE)
    plotname <- ''
    if(plot!=FALSE) {
      if(is.character(plot)) { plot <- paste0(plot,'.jpg'); jpeg(filename=plot); plotname <- plot }
      barplot(height=rbind(pj*100,kmj*100),beside=T,names.arg=1:length(pj),
              main=paste('Hosmer-Lemeshow (',version,') test\np-value =',format.pval(pval,digits=4),sep=''),
              xlab=xlab <- paste('Risk', tiles[length(pj)]),ylab='Risk, %',col=c('black','lightgray'))
      if(is.character(plot)) dev.off()
    }
    list(version=version,quantiles=tiles[length(pj)],nj=nj,kmj=kmj,pj=pj,chi2=dan_chi2,pval=pval,plotname=plotname)
  }
  h <- HosLem.test(surv_test, prediction)$pval
  
  # Add penalty in case of h <<< 0.05 (modify?)
  if (h < 1e-50) {
    c <- c - 0.3
  } else if (h > 1e-50 & h < 1e-10) {
    c <- c - 0.1
  } else if (h > 1e-10 & h < 0.05) {
    c <- c
  } else  if (h > 0.05) {
    c <- c + 0.05
  }
  
  return(c)
}

# Genetic algorithm
result <- ga(
  type = "binary",
  fitness = fitness,
  nBits = nbits,
  popSize = 100,
  pmutation = 0.3,
  pcrossover = 0.8,
  elitism = 5, 
  maxiter = 500,
  run = 50,
  parallel = 20,
  seed = 1993)


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


# Fit model
cvfit <- cv.glmnet(
  x = as.matrix(x_train),
  y = surv_train,
  type.measure = "C",
  family = "cox",
  alpha = 0)
prediction <- predict(cvfit,
                      as.matrix(subset(test, select = -c(Event_time, Event))),
                      s = "lambda.min")

# function to normalize prediction
normalize <- function(x, na.rm = TRUE) {
  return((x - min(x)) / (max(x) - min(x)))
}
prediction <- exp(prediction)
prediction <- as.vector(normalize(prediction))


# C-Index in test
c <- rcorr.cens(-prediction, surv_test)
print(paste0("C-Index in test set: ", c[1]))

# Ellapsed time
end <- Sys.time()
time <- difftime(end, start, units = "mins")
print(time)

# Save scores
# scores <- data.frame(
#   SampleID = rownames(test),
#   Score = prediction
# )
# print(head(scores))
# write.csv(scores, quote = FALSE, row.names = FALSE,
#           file = file.path(outputdir, "scores.csv"))