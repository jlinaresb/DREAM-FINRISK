get_scores <- function(coab_taxa, train, test, method) {

  taxa <- coab_taxa

  #TRAIN
  tt <- otu_table(train)@.Data
  tt <- tt[match(rownames(taxa), rownames(tt)), ]
  tt <- t(apply(tt, 1, function(x) log2(x + 1)))

  kk <- unique(taxa$clusters)
  signatures <- list()
  for (i in seq_along(kk)){
    signatures[[i]] <- rownames(taxa[which(taxa$clusters == kk[i]), ])
  }
  names(signatures) <- paste0("cluster_", kk)

  gsva_train <- gsva(tt,
                     signatures,
                     method = method,
                     kcdf = "Gaussian")
  gsva_train <- as.data.frame(t(gsva_train))

  #TEST
  tt <- otu_table(test)@.Data
  tt <- tt[match(rownames(taxa), rownames(tt)), ]
  tt <- t(apply(tt, 1, function(x) log2(x + 1)))

  gsva_test <- gsva(tt,
                    signatures,
                    method = method,
                    kcdf = "Gaussian")
  gsva_test <- as.data.frame(t(gsva_test))

  res <- list(train = gsva_train,
              test = gsva_test)

  return(res)
}