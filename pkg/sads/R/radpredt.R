radpredt <- function(object, ...){
  dots <- list(...)
  S <- length(object@data$x)
  if (object@sad == "ls" || object@sad == "geom" || object@sad == "nbinom"|| object@sad == "zipf"|| object@sad == "power"|| object@sad == "poilog"){
    y <- 1:max(object@data$x)
    if(!is.na(object@trunc)){
      if(object@sad=="ls") 
        X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(alpha = object@coef, N = sum(object@data$x)), lower.tail=F, trunc = object@trunc), dots))
      else
        X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
    }else {
      psad <- get(paste("p", object@sad, sep=""), mode = "function")
      if(object@sad=="ls")
        X <- do.call(psad, c(list(q = y, lower.tail = F, alpha = object@coef, N = sum(object@data$x)), dots))
      else
        X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(object@coef), dots))
    }
    f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
    ab <- f1(ppoints(S))
    if(is.na(ab[1]) & !any(is.na(ab[-1]))){
      ab[1] <- sum(object@data$x) - sum(ab[-1])
    }
  }else if(object@sad == "gamma" || object@sad == "lnorm" || object@sad == "weibull"){
    Y <- ppoints(S)
    if(!is.na(object@trunc)){
      ab <- do.call(qtrunc, c(list(object@sad, p = Y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
    }else{
      qsad <- get(paste("q", object@sad, sep=""), mode = "function")
      ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(object@coef), dots))
    }
  } else
    stop("unsupported distribution")
  new("rad", data.frame(rank=1:S, abund=ab))
}