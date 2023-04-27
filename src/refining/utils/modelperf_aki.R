# Copyright: Aki Havulinna @ THL

require(survival)
require(Hmisc) #cut2

# calculate absolute risk for a fitted Cox regression model, for "years" observation time
Coxar <- function(mod,years=10,bh=NULL) {
  if(is.null(bh)) bh <- basehaz(mod,centered=TRUE)
  stratum <- getStrata(mod,warn=FALSE)
  return(absrisk(surv=mod$y,bh=bh,lp=mod$linear.predictor,stratum=stratum,years=years))
}

# calculated absolute risk based coefficients of a fitted Cox regression model, but possibly newdata (the usual case)
Coxar.pred <- function(mod,betas=NULL,newdata=NULL,years=10) {
  if(is.null(betas)) betas <- coefficients(mod) # take the coefficients from the fitted model
  if(!is.null(newdata)) {
    xx <- model.matrix(mod,newdata)
    xx <- sweep(xx,2,colMeans(xx))
    if(!is.null(names(betas))) {
      if(!setequal(colnames(xx),names(betas))) stop('Coxar.pred - the names of model matrix and betas are different')
      betas <- betas[match(colnames(xx),names(betas))]
    }
    lp <- xx%*%betas
    stratum <- untangle.specials(mod$terms,'strata')$vars
    sterm <- paste(mod$formula[2]) # get the survival term
    fmula <- paste0(sterm,'~offset(lp)')
    if(length(stratum)>0) fmula <- paste(fmula,stratum,sep='+')
    mod <- coxph(as.formula(fmula),data=newdata)
  } else { # no newdata
    if(is.null(xx <- mod[['x']])) stop('either give data or use x=true in the model')
    xx <- sweep(xx,2,colMeans(xx))
    lp <- xx%*%betas
    stratum <- getStrata(mod,warn=FALSE)
    if(length(strat)>0) mod <- coxph(mod$y~strata(stratum)+offset(lp))
    else mod <- coxph(mod$y~offset(lp))
  }
  mod$offset <- NULL
  return(list(pred=c(Coxar(mod,years=years)),y=mod$y))
}

HosLem.test <- function(surv,pred,plot=FALSE,DF.reduce=2) {
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

#########################################
# the rest are necessary helper functions
#########################################

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


