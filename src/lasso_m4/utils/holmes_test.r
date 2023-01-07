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