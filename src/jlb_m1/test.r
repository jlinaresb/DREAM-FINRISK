# Validation
# ========
setwd(here::here())
source("src/utils/importPseq.r")
source("src/jlb_m1/preprocessing.r")
source("requirements.r")

models <- readRDS("src/jlb_m1/results/models.rds")

# 
models$train <- models$train[complete.cases(models$train), ]

test <- pseq(subset = "test")
test <- remove_samples(test)
test <- norm(test, filter = FALSE)

taxids <- colnames(models$train)
test <- prune_taxa(taxids, test)

otu_test <- as.data.frame(t(otu_table(test)@.Data))
pheno_test <- sample_data(test)
data_test <- cbind.data.frame(pheno_test, otu_test)

data_test <- data_test[which(!is.na(data_test$Event)), ]
data_impute <- impute.knn(t(data_test))
data_test <- t(data_impute$data)

predRes(res = models$fit,
        method = "ridgelasso",
        traindata = models$train,
        newdata = data_test,
        int.cv = FALSE,
        int.cv.nfold = 5,
        time = seq(2, 15, 1),
        trace = TRUE,
        ncores = 8
        )
