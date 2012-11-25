qqsad <- function(object, sad, coef, trunc=NA, distr){
  if(class(object)=="fitsad"){
    sad <- object@sad
    coef <- as.list(bbmle::coef(object))
    trunc <- object@trunc
    distr <- object@distr
    x <- object@data$x
  }
  else if(class(object)=="numeric")
    x <- object
  rank <- sort(x)
  S <- length(x)
  if(distr == "D"){
    q <- 1:sum(x)
    if(!is.na(trunc)){
      if(sad == "ls")
        p <- do.call(ptrunc,
                     c(list(sad, q = q, coef = c(sum(x), as.numeric(coef)), trunc = trunc)))
      else if(sad == "volkov")
        p <- do.call(ptrunc,
                     c(list(sad, q = q, coef = c(as.numeric(coef), sum(x)), trunc = trunc)))
      else if(sad == "mzsm")
        p <- do.call(ptrunc,
                     c(list(sad, q = q, coef = c(as.numeric(coef), sum(x)), trunc = trunc)))
      else
        p <- do.call(ptrunc, list(sad, q = q, trunc = trunc, coef=coef))
    }
    else{
      if(sad == "ls")
        p <- do.call(pls, c(list(q = q), N = sum(x), alpha = as.numeric(coef)))
      else if(sad =="volkov")
        p <- do.call(pvolkov, c(list(q = q, J=sum(x)), coef))
      else if(sad =="mzsm")
        p <- do.call(pmzsm, c(list(q = q, J=sum(x)), coef))
      else{
        psad <- get(paste("p", sad, sep=""), mode = "function")
        p <- do.call(psad, c(list(q = q), coef))
      }
    }
    f1 <- approxfun(x=c(1, p), y=c(0, q), method="constant")
    q <- f1(ppoints(S))
  }
  else if(distr == "C"){
    p <- ppoints(S)
    if(!is.na(trunc))
      q <- do.call(qtrunc, list(sad, p = p, trunc = trunc,coef=coef))
    else{
      qsad <- get(paste("q", sad, sep=""), mode = "function")
      q <- do.call(qsad, c(list(p = p), coef))
    }
  }
  else
    stop("unsupported distribution")
  plot(q, rank, main = "Q-Q plot", xlab="Theoretical Quantile", ylab="Sample Quantiles")
  abline(0, 1, col = "red", lty = 2)
}
