
# (C) Aki Havulinna, 2022-06-17
#library(Hmisc) # if it works, otherwisae just use cut2 given below

# Coxar and Hoslem.test are the functions of interest
# the rest are helper functions 

tiles <- c('','','tertiles','quartiles','quintiles','sextiles','septiles','octiles','noniles','deciles')

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

 # Below, pred is the predicted absolute risk (e.g. by Coxar)
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
  grp <- as.numeric(cut2(pred, g=10))
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


getStrata <- function(coxfit,newdata=NULL,warn=TRUE) {
  # this function extracts strata for each observation in the coxfit (or model frame) object
  # safer to use "newdata"
  temp <- untangle.specials(terms(coxfit), "strata")
  if(length(temp$vars)<1) {
    if(warn) warning('no strata in the model')
    return(NULL)
  }
  if(!is.null(newdata)) mf <- model.frame(coxfit,data=newdata)
  else mf <- model.frame(coxfit)
  if (length(temp$vars) == 1) stratum <- mf[[temp$vars]]
  else stratum <- strata(mf[, temp$vars], shortlabel = TRUE)
  stratum
}

findhaz <- function(bh,time0,time1=NULL,stratum=NULL) {
  if(is.null(bh$strata)) {
    if(!is.null(stratum)) stop('did not expect strata')
    af <- approxfun(bh$time,bh$hazard,rule=2)
    h <- -af(time0)
    if(!is.null(time1)) h <- -h-af(time1)
  } else {
    if(is.null(stratum)) stop('basehaz has strata, expected strata as an argument')
    h <- rep(NA,length(time0))
    us <- unique(bh$strata)
    if(length(setdiff(us,unique(stratum)))) stop('stratum has categories which are not in baseline')
    for(st in us) {
      idx <- bh$strata==st; idx2 <- stratum==st
      af <- approxfun(bh$time[idx],bh$hazard[idx],rule=2)
      h[idx2] <- -af(time0[idx2]) # this is h(end) if time1 is not given
      if(!is.null(time1)) h[idx2] <- -h[idx2]-af(time1[idx2])
    }
  }
  h
}

absrisk <- function(surv,bh,lp=0,stratum=NULL,years=10) {
  # calculates absolute risk given survival indicators, baseline hazard and linear predictor + possible strata
  n <- nrow(surv)
  type <- attr(surv,'type')
  if(type=='counting') {
    haz <- findhaz(bh,surv[,1],surv[,1]+years,stratum=stratum)
  } else if (type=='right') {
    haz <- findhaz(bh,rep(years,n),stratum=stratum)
  } else {
    stop('Unsupported censoring type: ',type)
  }
  1-exp(haz*exp(lp))
}

# pedicted absolute risk, from Cox model
Coxar <- function(mod,years,bh=NULL) {
  if(is.null(bh)) bh <- basehaz(mod,centered=TRUE)
  stratum <- getStrata(mod,warn=FALSE)
  return(absrisk(surv=mod$y,bh=bh,lp=mod$linear.predictor,stratum=stratum,years=years))
}


# cut2 is borrowed from Hmisc, which does not seem to install on Atlas.
cut2 <- function (x, cuts, m = 150, g, levels.mean = FALSE, digits, minmax = TRUE, 
          oneval = TRUE, onlycuts = FALSE, formatfun = format, ...) 
{
  if (inherits(formatfun, "formula")) {
    if (!requireNamespace("rlang")) 
      stop("Package 'rlang' must be installed to use formula notation")
    formatfun <- getFromNamespace("as_function", "rlang")(formatfun)
  }
  method <- 1
  x.unique <- sort(unique(c(x[!is.na(x)], if (!missing(cuts)) cuts)))
  min.dif <- min(diff(x.unique))/2
  min.dif.factor <- 1
  if (missing(digits)) 
    digits <- if (levels.mean) 
      5
  else 3
  format.args <- if (any(c("...", "digits") %in% names(formals(args(formatfun))))) {
    c(digits = digits, list(...))
  }
  else {
    list(...)
  }
  oldopt <- options("digits")
  options(digits = digits)
  on.exit(options(oldopt))
  xlab <- attr(x, "label")
  if (missing(cuts)) {
    nnm <- sum(!is.na(x))
    if (missing(g)) 
      g <- max(1, floor(nnm/m))
    if (g < 1) 
      stop("g must be >=1, m must be positive")
    options(digits = 15)
    n <- table(x)
    xx <- as.double(names(n))
    options(digits = digits)
    cum <- cumsum(n)
    m <- length(xx)
    y <- as.integer(ifelse(is.na(x), NA, 1))
    labs <- character(g)
    cuts <- approx(cum, xx, xout = (1:g) * nnm/g, method = "constant", 
                   rule = 2, f = 1)$y
    cuts[length(cuts)] <- max(xx)
    lower <- xx[1]
    upper <- 1e+45
    up <- low <- double(g)
    i <- 0
    for (j in 1:g) {
      cj <- if (method == 1 || j == 1) 
        cuts[j]
      else {
        if (i == 0) 
          stop("program logic error")
        s <- if (is.na(lower)) 
          FALSE
        else xx >= lower
        cum.used <- if (all(s)) 
          0
        else max(cum[!s])
        if (j == m) 
          max(xx)
        else if (sum(s) < 2) 
          max(xx)
        else approx(cum[s] - cum.used, xx[s], xout = (nnm - 
                                                        cum.used)/(g - j + 1), method = "constant", 
                    rule = 2, f = 1)$y
      }
      if (cj == upper) 
        next
      i <- i + 1
      upper <- cj
      y[x >= (lower - min.dif.factor * min.dif)] <- i
      low[i] <- lower
      lower <- if (j == g) 
        upper
      else min(xx[xx > upper])
      if (is.na(lower)) 
        lower <- upper
      up[i] <- lower
    }
    low <- low[1:i]
    up <- up[1:i]
    variation <- logical(i)
    for (ii in 1:i) {
      r <- range(x[y == ii], na.rm = TRUE)
      variation[ii] <- diff(r) > 0
    }
    if (onlycuts) 
      return(unique(c(low, max(xx))))
    flow <- do.call(formatfun, c(list(low), format.args))
    fup <- do.call(formatfun, c(list(up), format.args))
    bb <- c(rep(")", i - 1), "]")
    labs <- ifelse(low == up | (oneval & !variation), flow, 
                   paste("[", flow, ",", fup, bb, sep = ""))
    ss <- y == 0 & !is.na(y)
    if (any(ss)) 
      stop(paste("categorization error in cut2.  Values of x not appearing in any interval:\n", 
                 paste(format(x[ss], digits = 12), collapse = " "), 
                 "\nLower endpoints:", paste(format(low, digits = 12), 
                                             collapse = " "), "\nUpper endpoints:", paste(format(up, 
                                                                                                 digits = 12), collapse = " ")))
    y <- structure(y, class = "factor", levels = labs)
  }
  else {
    if (minmax) {
      r <- range(x, na.rm = TRUE)
      if (r[1] < cuts[1]) 
        cuts <- c(r[1], cuts)
      if (r[2] > max(cuts)) 
        cuts <- c(cuts, r[2])
    }
    l <- length(cuts)
    k2 <- cuts - min.dif
    k2[l] <- cuts[l]
    y <- cut(x, k2)
    if (!levels.mean) {
      brack <- rep(")", l - 1)
      brack[l - 1] <- "]"
      fmt <- do.call(formatfun, c(list(cuts), format.args))
      labs <- paste("[", fmt[1:(l - 1)], ",", fmt[2:l], 
                    brack, sep = "")
      if (oneval) {
        nu <- table(cut(x.unique, k2))
        if (length(nu) != length(levels(y))) 
          stop("program logic error")
        levels(y) <- ifelse(nu == 1, c(fmt[1:(l - 2)], 
                                       fmt[l]), labs)
      }
      else levels(y) <- labs
    }
  }
  if (levels.mean) {
    means <- tapply(x, y, function(w) mean(w, na.rm = TRUE))
    levels(y) <- do.call(formatfun, c(list(means), format.args))
  }
  attr(y, "class") <- "factor"
  if (length(xlab)) 
    label(y) <- xlab
  y
}
