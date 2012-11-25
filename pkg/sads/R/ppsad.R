ppsad <- function (object, sad, coef, trunc=NA, distr) {
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
  z <- ppoints(S)
  if(!is.na(trunc)){
    if(sad == "ls")
      p <- do.call(ptrunc, c(list(sad, q = rank, coef = c(sum(x), as.numeric(coef)), trunc = trunc)))
    else if(sad == "volkov")
      p <- do.call(ptrunc, c(list(sad, q = rank, coef = c(as.numeric(coef),sum(x))), trunc = trunc))
    else if(sad == "mzsm")
      p <- do.call(ptrunc, c(list(sad, q = rank, coef = c(as.numeric(coef),sum(x))), trunc = trunc))
    else
      p <- do.call(ptrunc, list(sad, q = rank, coef = coef, trunc = trunc))
  }
  else{
    if(sad == "ls")
      p <- do.call(pls, c(list(q = rank), N=sum(x), alpha=as.numeric(coef)))
    else if(sad == "volkov")
      p <- do.call(pvolkov, c(list(q = rank, J=sum(x)), as.numeric(coef)))
    else if(sad == "mzsm")
      p <- do.call(pmzsm, c(list(q = rank, J=sum(x)), as.numeric(coef)))
    else{
      psad <- get(paste("p", sad, sep=""), mode = "function")
      p <- do.call(psad, c(list(q = rank), coef))
    }
  }
  plot(z, p, main = "P-P plot", ylim = c(0, 1), xlab='Theoretical Percentiles', ylab='Sample Percentiles')
  abline(0, 1, col = "red", lty = 2)
}
