dataTrans <- function(data, x, y, z, tt, std.x, std.i, std.tt, inter, trace = TRUE){
  
  #########################################################
  ### Control and checking
  
  if(length(y) != 2)
    stop("\nThe outcome must contains 2 columns: 'time' and 'status'.")
  if(min(data[, y[1]]) < 0)
    stop("\nTimes must be non-negative.")
  if(sum(data[, y[2]] %in% c(0, 1)) != nrow(data))
    stop("\nStatus must be either 0 (censor) or 1 (event).")
  
  if(!is.null(tt)){
    TT <- as.data.frame(data[, tt])
    if(length(tt) > 1)
      stop("\nThe treatment must be a single variable.")
    if(length(tt) == 1){
      if(length(unique(as.vector(t(data[, tt])))) > 2)
        stop(paste0("\nThe treatment variable must consider only two groups."))
      if(length(unique(as.vector(t(data[, tt])))) == 1){
        warning(paste0("\nAll patients are in the same treatment group. The analysis is then switch to a prognostic setting."))
        tt <- NULL
        inter <- FALSE
      }
      if((sum(unique(as.vector(t(data[, tt]))) %in% c(-0.5, +0.5)) != 2) & std.tt == TRUE)
        data[, tt] <- as.numeric(factor(as.vector(t(data[, tt])))) - 1.5
    }
  }
  
  #########################################################
  
  itt <- tt; iz <- z; ix <- x; iy <- y
  iptt <- which(colnames(data) %in% tt)
  ipz <- which(colnames(data) %in% z)
  ipx <- which(colnames(data) %in% x)
  ipy <- which(colnames(data) %in% y)
  isSim <- (!is.null(attributes(data)$isSim))
  
  if(std.x == TRUE)
    data[, x] <- scale(data[, x], center = T, scale = T)
  
  if(inter == TRUE){
    XT <- as.matrix(data[, x]) * matrix(data[, tt], nrow = nrow(data), ncol = length(x))
    if(std.i == TRUE)
      XT <- scale(XT, center = T, scale = T)
    colnames(XT) <- paste0("bi", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
    xt <- colnames(XT)
  }else{
    xt <- NULL
  }
  
  data <- cbind(data[, c(tt, z, x, y)])
  if(inter == TRUE){
    data <- cbind(data, XT)
    colnames(data)[1] <- "treat"
  }
  
  tnames <- c(rep("tt", length(tt)), rep("z", length(z)), rep("x", length(x)), rep("y", length(y)), rep("xt", length(xt)))
  if(!is.null(z))
    colnames(data)[tnames == "z"] <- paste0("cl", gsub(" ", "0", format(c(length(z), 1:length(z))))[-1])
  colnames(data)[tnames == "x"] <- paste0("bm", gsub(" ", "0", format(c(length(x), 1:length(x))))[-1])
  colnames(data)[tnames == "y"] <- c("time", "status")
  
  data <- na.omit(data)
  if(!is.null(attributes(data)$na.action) & trace == TRUE){
    nmiss <- length(attributes(data)$na.action)
    message(paste0(
      "\rData management: ", nmiss, " observation", ifelse(nmiss > 1, "s were", " was"), " excluded due to missing data."))
  }
  
  if(!(class(unlist(data)) %in% c("numeric", "integer")))
    stop("\nAll variables must be numerical.")
  
  attributes(data) <- append(
    x = attributes(data),
    values = list(
      inter = inter,
      inames = list(
        tt = itt,
        z = iz,
        x = ix,
        y = iy),
      ipos = list(
        tt = iptt,
        z = ipz,
        x = ipx,
        y = ipy
      ),
      tnames = tnames,
      pos = list(
        z = grep("cl", colnames(data)),
        x = grep("bm", colnames(data)),
        xt = grep("bi", colnames(data)),
        X = (1:ncol(data))[-which(colnames(data) %in% c("time", "status"))],
        y = which(colnames(data) %in% c("time", "status"))
      ),
      weights = c(rep(0, length(c(tt, z))), rep(1, length(x) * ((inter == TRUE) + 1))),
      std.x = std.x,
      std.tt = std.tt)
  )
  if(inter == TRUE){
    attributes(data)$pos$tt <- grep("treat", colnames(data))
    attributes(data)$std.i <- std.i
    attributes(data)$inames$xt <- paste0(ix, ":", itt)
  }
  
  return(data)
}

predRes2 <- function(
    ####################################################################
    ######################## *** PARAMETERS *** ########################
    res,              # Object of class 'resBMsel'
    method,           # Methods to compute
    traindata,        # Training dataset
    newdata,          # New dataset
    int.cv,           # Internal CV should be performed?
    int.cv.nfold = 5, # Number of folds for the internal CV
    time,             # Time points
    trace = TRUE,     # Print function's progression?
    ncores = 1        # Number of PC cores used
    ####################################################################
){
  
  ####################################################################
  ###  DATA CHECKING AND MANIPULATION
  
  if(class(res) != "resBMsel")
    stop("\n'res' must be an object returned by the function BMsel().")
  
  if(missing(method)){
    method <- colnames(summary(res, show = FALSE, add.ridge = !is.na(attributes(res)$ridge)))
  }else{
    if(length(setdiff(method, names(res))) > 0)
      stop("\n Some methods in 'method' were not previously computed or do not exist.")
    method <- unique(c(method, if(attributes(res)$isSim == TRUE) "oracle"))
  }
  
  if(missing(traindata))
    stop("\nThe training data set used in the BMsel() function must be specified.")
  traindata <- as.data.frame(traindata)
  rownames(traindata) <- 1:nrow(traindata)
  
  if(length(setdiff(attributes(res)$inames[which(attributes(res)$tnames != 'xt')], colnames(traindata))) > 0)
    stop("\nSome covariates of the training set are missing. Please specify the same data set used for BMsel().")
  
  if(missing(int.cv))
    int.cv <- FALSE
  
  if(!(int.cv) %in% c(TRUE, FALSE))
    stop("\n'int.cv' must be either TRUE or FALSE.")
  
  if(int.cv.nfold < 2 || int.cv.nfold > nrow(traindata))
    stop("\n'int.cv.nfold' must be between 2 and the sample size 'n'.")
  
  # if(missing(time)){
  #   stop("\n'time' must be specified for Cox models.")
  # }else{
  #   if(min(time) < 0 || max(time) > min(c(max(traindata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')[1]]]),
  #                                         if(!missing(newdata)) max(newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')[1]]]))))
  #     stop("\n'time' is out of the range of the observed survival time.")
  # }
  time <- sort(time)
  
  # if(!missing(newdata)){
  #   if(length(setdiff(attributes(res)$inames[which(attributes(res)$tnames != 'xt')], colnames(newdata))) > 0)
  #     stop("Some covariates of the new data set are missing. Please specify the same covariates as the training set.")
  #   newdata <- as.data.frame(newdata)
  #   rownames(newdata) <- 1:nrow(newdata)
  # }
  
  ncores <- round(ncores, 0)
  if(ncores < 1 || ncores > parallel::detectCores())
    stop(paste0("\n'ncores' must be between 1 and ", parallel::detectCores(), "."))
  if(ncores > int.cv.nfold) ncores <- int.cv.nfold
  
  tt <- attributes(res)$inames[which(attributes(res)$tnames == 'tt')]
  x <- attributes(res)$inames[which(attributes(res)$tnames == 'x')]
  z <- attributes(res)$inames[which(attributes(res)$tnames == 'z')]
  y <- attributes(res)$inames[which(attributes(res)$tnames == 'y')]
  isRidge <- (unique(!is.na(attributes(res)$ridge)))
  isNew <- (!missing(newdata))
  
  Res <- data.frame(summary(res, show = FALSE, add.ridge = isRidge))
  Res <- Res[, gsub("-", ".", method), drop = FALSE]
  
  if(attributes(res)$inter == TRUE){
    Res.i <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
    Res.i <- Res.i[, gsub("-", ".", method), drop = FALSE]
  }
  
  nmeth <- ncol(Res)
  
  ####################################################################
  ### PREDICTION FOR TRAINING SET
  
  # if(trace == TRUE)
  #   message(paste0(
  #     "\rComputing prediction criteria for: training set"))
  # 
  # attr.train <- attributes(traindata)[-which(names(attributes(traindata)) %in% "names")]
  # 
  # tdata <- dataTrans(data = traindata, x = x, y = y, z = z, tt = tt,
  #                    std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
  #                    inter = attributes(res)$inter, trace = FALSE)
  # colnames(tdata) <- attributes(res)$inames
  # attributes(tdata) <- append(attributes(tdata), attr.train)
  # 
  # surv.train <- Surv(tdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][1]],
  #                    tdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][2]])
  # 
  # lp.train <- matrix(0, nrow = nrow(tdata), ncol = ncol(Res))
  # if(nrow(Res) > 0)
  #   lp.train <- as.matrix(tdata[, rownames(Res)]) %*% as.matrix(Res)
  # 
  # if(attributes(res)$inter == TRUE){
  #   lpint.train <- matrix(0, nrow = nrow(tdata), ncol = ncol(Res))
  #   if(nrow(Res.i) > 0 & sum(Res.i) != 0)
  #     lpint.train <- as.matrix(
  #       tdata[, gsub(paste0(":", attributes(res)$inames[1]), "", rownames(Res.i))]) %*% as.matrix(Res.i)
  # }
  # 
  # predRes.train <- compute.predRes(res = res, nmeth = nmeth, hrz = time, traindata = traindata, newdata = traindata,
  #                                  surv.train = surv.train, surv.new = surv.train, lp.train = lp.train, lp.new = lp.train,
  #                                  lpint.train = lpint.train, lpint.new = lpint.train, tt = tt)
  
  
  ####################################################################
  ### PREDICTION FOR VALIDATION SET
  
  if(!missing(newdata)){
    if(trace == TRUE)
      message(paste0(
        "\rComputing prediction criteria for: validation set"))
    
    newdata <- dataTrans(data = newdata, x = x, y = y, z = z, tt = tt,
                         std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
                         inter = attributes(res)$inter, trace = TRUE)
    colnames(newdata) <- attributes(res)$inames
    
    # surv.new <- Surv(newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][1]],
    #                  newdata[, attributes(res)$inames[which(attributes(res)$tnames == 'y')][2]])
    
    lp.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
    if(nrow(Res) > 0)
      lp.new <- as.matrix(newdata[, rownames(Res)]) %*% as.matrix(Res)
    
    if(attributes(res)$inter == TRUE){
      lpint.new <- matrix(0, nrow = nrow(newdata), ncol = ncol(Res))
      if(nrow(Res.i) > 0 & sum(Res.i) != 0)
        lpint.new <- as.matrix(
          newdata[, gsub(paste0(":", attributes(res)$inames[1]), "", rownames(Res.i))]) %*% as.matrix(Res.i)
    }
    
    # predRes.new <- compute.predRes(res = res, nmeth = nmeth, hrz = time, traindata = traindata, newdata = newdata,
    #                                surv.train = surv.train, surv.new = surv.new, lp.train = lp.train, lp.new = lp.new,
    #                                lpint.train = lpint.train, lpint.new = lpint.new, tt = tt)
  }
  
  
  ####################################################################
  ### PREDICTION FOR INTERNAL DOUBLE CROSS-VALIDATION (2CV)
  
  # if(int.cv == TRUE){
  #   form <- attributes(res)$formula
  #   foldid2 <- sample(x = 1:int.cv.nfold, size = nrow(traindata), replace = T)
  #   
  #   if(trace == TRUE)
  #     message(
  #       "\rComputing prediction criteria for: internal validation")
  #   
  #   cl <- makeCluster(ncores)
  #   
  #   res.int.cv <- clusterApplyLB(
  #     cl = cl,
  #     x = 1:int.cv.nfold,
  #     fun = function(X){
  #       
  #       traindataT <- traindata[which(foldid2 != X), ]
  #       ltraindataT <- list(traindataT)
  #       traindataV <- traindata[which(foldid2 == X), ]
  #       
  #       form[2] <- ltraindataT
  #       form[which(names(form) == "method")] <- list(setdiff(method, "oracle"))
  #       w <- options()$warn
  #       options(warn = -1)
  #       pos.trace <- which(names(form) == "trace")
  #       if(length(pos.trace) == 0){
  #         form[length(names(form)) + 1] <- FALSE
  #         names(form)[length(names(form))] <- "trace"
  #       }else{
  #         form[pos.trace] <- FALSE
  #       }
  #       res <- eval(form)
  #       options(warn = w)
  #       
  #       Res <- data.frame(summary(res, show = FALSE, add.ridge = isRidge))
  #       Res <- Res[, gsub("-", ".", method), drop = FALSE]
  #       
  #       if(attributes(res)$inter == TRUE){
  #         Res.i <- data.frame(summary(res, show = FALSE, keep = "xt", add.ridge = isRidge))
  #         Res.i <- Res.i[, gsub("-", ".", method), drop = FALSE]
  #       }
  #       
  #       nmeth <- ncol(Res)
  #       
  #       traindataV <- dataTrans(data = traindataV, x = x, y = y, z = z, tt = tt,
  #                               std.x = attributes(res)$std.x, std.i = attributes(res)$std.i, std.tt = attributes(res)$std.tt,
  #                               inter = attributes(res)$inter, trace = FALSE)
  #       colnames(traindataV) <- attributes(res)$inames
  #       
  #       lp.trainV <- matrix(0, nrow = nrow(traindataV), ncol = ncol(Res))
  #       if(nrow(Res) > 0)
  #         lp.trainV <- as.matrix(traindataV[, rownames(Res)]) %*% as.matrix(Res)
  #       
  #       lpint.trainV <- NA
  #       if(attributes(res)$inter == TRUE){
  #         lpint.trainV <- matrix(0, nrow = nrow(traindataV), ncol = ncol(Res))
  #         if(nrow(Res.i) > 0 & sum(Res.i) != 0)
  #           lpint.trainV <- as.matrix(
  #             traindataV[, gsub(paste0(":", attributes(res)$inames[1]), "", rownames(Res.i))]) %*% as.matrix(Res.i)
  #       }
  #       
  #       rownames(lp.trainV) <- rownames(traindataV)
  #       colnames(lp.trainV) <- method
  #       
  #       if(attributes(res)$inter == TRUE){
  #         rownames(lpint.trainV) <- rownames(traindataV)
  #         colnames(lpint.trainV) <- method
  #       }
  #       
  #       return(list(lp.trainV, lpint.trainV))
  #     }
  #   )
  #   
  #   stopCluster(cl)
  #   
  #   lp.int.cv <- data.frame()
  #   for(i in 1:int.cv.nfold)
  #     lp.int.cv <- rbind(lp.int.cv, data.frame(res.int.cv[[i]][1]))
  #   lp.int.cv <- lp.int.cv[order(as.numeric(rownames(lp.int.cv))), , drop = FALSE]
  #   
  #   if(attributes(res)$inter == TRUE){
  #     lpint.int.cv <- data.frame()
  #     for(i in 1:int.cv.nfold)
  #       lpint.int.cv <- rbind(lpint.int.cv, data.frame(res.int.cv[[i]][2]))
  #     lpint.int.cv <- lpint.int.cv[order(as.numeric(rownames(lpint.int.cv))), , drop = FALSE]
  #   }
  #   
  #   predRes.int.cv <- compute.predRes(res = res, nmeth = nmeth, hrz = time, traindata = traindata, newdata = traindata,
  #                                     surv.train = surv.train, surv.new = surv.train, lp.train = lp.train, lp.new = lp.int.cv,
  #                                     lpint.train = lpint.train, lpint.new = lpint.int.cv, tt = tt)
  #   
  # }
  
  
  ####################################################################
  ### FORMATTING RESULTS
  
  # predRes <- list()
  # for(i in 1:length(time)){
  #   predres <- list('Training set' = round(predRes.train[[i]], 4))
  #   if(int.cv == TRUE) predres <- merge.list(predres, list('Internal validation' = round(predRes.int.cv[[i]], 4)))
  #   if(!missing(newdata)) predres <- merge.list(predres, list('External validation' = round(predRes.new[[i]], 4)))
  #   predRes[[i]] <- predres
  # }
  # names(predRes) <- paste0("time = ", time)
  # class(predRes) <- "predRes"
  return(lp.new)
}