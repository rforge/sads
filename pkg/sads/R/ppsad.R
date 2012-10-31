ppsad <- function (object) {
  rank <- sort(object@data$x)
  S <- length(object@data$x)
  z <- ppoints(S)
  if(!is.na(object@trunc)){
    if(object@sad == "ls")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
    else if(object@sad == "zipf")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(length(object@data$x), as.numeric(object@coef))), trunc = object@trunc))
    else if(object@sad == "volkov")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(sum(object@data$x), as.numeric(object@coef))), trunc = object@trunc))
    else if(object@sad == "mzsm")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(sum(object@data$x), as.numeric(object@coef))), trunc = object@trunc))
    else if(object@sad == "mand")
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = c(sum(object@data$x), as.numeric(object@coef))), trunc = object@trunc))
    else
      p <- do.call(ptrunc, c(list(object@sad, q = rank, coef = as.list(object@coef), trunc = object@trunc)))
  }else{
    if(object@sad == "ls")
      p <- do.call(pls, c(list(q = rank), sum(object@data$x), as.numeric(object@coef)))
    else if(object@sad == "zipf")
      p <- do.call(pzipf, c(list(q = rank), length(object@data$x), as.numeric(object@coef)))
    else if(object@sad == "volkov")
      p <- do.call(pvolkov, c(list(q = rank), sum(object@data$x), as.numeric(object@coef)))
    else if(object@sad == "mzsm")
      p <- do.call(pmzsm, c(list(q = rank), sum(object@data$x), as.numeric(object@coef)))
    else if(object@sad == "mand")
      p <- do.call(pmand, c(list(q = rank), sum(object@data$x), as.numeric(object@coef)))
    else{
      psad <- get(paste("p", object@sad, sep=""), mode = "function")
      p <- do.call(psad, c(list(q = rank), as.list(object@coef)))
    }
  }
  plot(z, p, main = "P-P plot", ylim = c(0, 1), xlab='Theoretical Percentiles', ylab='Sample Percentiles')
  abline(0, 1, col = "red", lty = 2)
}