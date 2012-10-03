plotprofmle <- function(mleobj, nseg=20, ratio=log(8), which=NULL, auto.mfrow=TRUE,
                        ask=!auto.mfrow&(prod(par("mfcol")) < length(which) && dev.interactive()),
                        col.line="blue", varname, ...){
  if( class(mleobj)[1] != "profile.mle" &
     class(mleobj)[1] != "profile.mle2") 
    stop( "Object should have class \'profile.mle\' or \'profile.mle2\'")
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  mleprof <- mleobj@profile
  npar <- length(mleprof)
  dots <- list(...)
  if(!"ylab" %in% names(dots)) dots$ylab <- "Negative relative log-likelihood"
  if(!"type" %in% names(dots)) dots$type <- "l"
  if(!"col" %in% names(dots)) dots$col <- "red"
  vname <- names(mleprof)
  if( is.null(which) ){
    if(!missing(varname)){
      if(length(varname)!=length(mleprof))stop("Length of 'varname' should match number os mles in mle.prof object")
      vname <- varname
    }
    parseq = 1:npar
    if (auto.mfrow) {
      nl <- floor(sqrt(npar))
      nc <- ceiling(npar/nl)
      par(mfrow=c(nl,nc))
    }
  }
  else{
    if(!missing(varname)){
      if(length(varname)!=length(which))stop("Length of 'which' should match length of 'varnames'")
      vname[which] <- varname
    }
    parseq = which
  }
  if (ask&!auto.mfrow) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  for(i in parseq)
    {
      tmp <- mleprof[i][[1]]
      y <- tmp[,1]^2/2
      x <- (tmp[,2][,i])
      interpol = spline(x, y, n=nseg*length(x) )
      do.call(plot, c(list(x=interpol,xlab=vname[i]),dots))
      if(!is.null(ratio)){
        l <- length(interpol$y)
        change <- (interpol$y - ratio)[2:l] * (interpol$y - ratio)[1:(l-1)]
        endpoints <- which(change < 0)
        corr <- (interpol$x[2]-interpol$x[1])/2
        for (j in 1:(length(endpoints)/2)) {
          lower <-interpol$x[endpoints[(2*j)-1]]+corr
          upper <- interpol$x[endpoints[2*j]]+corr
          lines(c(lower,upper ),c(ratio, ratio), col=col.line, lty=2)
          lines(rep(lower,2), c(-1, ratio), col=col.line, lty=2)
          lines(rep(upper,2), c(-1, ratio), col=col.line, lty=2)
        }
      }
    }
}
