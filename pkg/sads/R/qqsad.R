qqsad <- function(object){
  rank <- sort(object@data$x)
  S <- length(object@data$x)
  if(object@distr == "D"){
    q <- 1:sum(object@data$x)
    if(!is.na(object@trunc)){
      if(object@sad == "ls")
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
      else if(object@sad == "zipf")
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = c(length(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
      else if(object@sad == "volkov")
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
      else if(object@sad == "mzsm")
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
      else if(object@sad == "mand")
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = c(sum(object@data$x), as.numeric(object@coef)), trunc = object@trunc)))
      else
        p <- do.call(ptrunc, c(list(object@sad, q = q, coef = as.list(object@coef), trunc = object@trunc)))
    }else{
      if(object@sad == "ls")
        p <- do.call(pls, c(list(q = q), N = sum(object@data$x), alpha = as.numeric(object@coef)))
      else if(object@sad == "zipf")
        p <- do.call(pzipf, c(list(q = q), N = length(object@data$x), s = as.numeric(object@coef)))
      else if(object@sad =="volkov")
        p <- do.call(pvolkov, c(list(q = q), sum(object@data$x), as.numeric(object@coef)))
      else if(object@sad =="mzsm")
        p <- do.call(pmzsm, c(list(q = q), sum(object@data$x), as.numeric(object@coef)))
      else if(object@sad =="mand")
        p <- do.call(pmand, c(list(q = q), sum(object@data$x), as.numeric(object@coef)))
      else{
        psad <- get(paste("p", object@sad, sep=""), mode = "function")
        p <- do.call(psad, c(list(q = q), as.list(object@coef)))
      }
    }
    f1 <- approxfun(x=c(1, p), y=c(0, q), method="constant")
    q <- f1(ppoints(S))
  }else if(object@distr == "C"){
    p <- ppoints(S)
    if(!is.na(object@trunc))
      q <- do.call(qtrunc, c(list(object@sad, p = p, coef = as.list(object@coef), trunc = object@trunc)))
    else{
      qsad <- get(paste("q", object@sad, sep=""), mode = "function")
      if(object@sad == "pareto")
        q <- do.call(qsad, c(list(p = p), shape = as.numeric(object@coef), scale=min(object@data$x)))
      else
        q <- do.call(qsad, c(list(p = p), as.list(object@coef)))
    }
  }else
    stop("unsupported distribution")
  plot(q, rank, main = "Q-Q plot", xlab="Theoretical Quantile", ylab="Sample Quantiles")
  abline(0, 1, col = "red", lty = 2)
}