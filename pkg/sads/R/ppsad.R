ppsad <- function (object) {
  rank <- sort(object@data$x)
  S <- length(object@data$x)
  z <- ppoints(S)
  if(!is.na(object@trunc)){
    if(object@sad == "ls")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
    else
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = as.list(object@coef), trunc = object@trunc)))
  }else{
    if(object@sad == "ls")
      p <- do.call(pls, c(list(q = rank), sum(object@data$x), as.numeric(object@coef)))
    else{
      psad <- get(paste("p", object@sad, sep=""), mode = "function")
      p <- do.call(psad, c(list(q = rank), as.list(object@coef)))
    }
  }
  plot(z, p, main = "P-P plot", ylim = c(0, 1), xlab='Theoretical Percentiles', ylab='Sample Percentiles')
  abline(0, 1, col = "red", lty = 2)
}