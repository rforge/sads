radpred <- function(object, x, sad, coef, trunc, S, A, ...){
  dots <- list(...)
  if(!missing(sad) && !missing(coef)){
    if (!missing(x)){
      S <- length(x)
      if (sad=="mand" || sad == "volkov" || sad == "mzsm" || sad == "ls" || sad == "geom" || sad == "nbinom"|| sad == "zipf"|| sad == "power"|| sad == "poilog"){
        y <- 1:max(x)
        if(!missing(trunc)){
          X <- do.call(ptrunc, c(list(sad, q = y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
        }else{
          psad <- get(paste("p", sad, sep=""),  mode = "function")
          X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(coef), dots))
        }
        f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
        ab <- f1(ppoints(S))
        #if(!any(is.na(ab[-1]))){
        #  ab[1] <- sum(x) - sum(ab[-1])
        #}
      }else if(sad == "pareto" || sad == "gamma" || sad == "lnorm" || sad == "weibull"){
        Y <- ppoints(S)
        if(!missing(trunc)){
          ab <- do.call(qtrunc, c(list(sad, p = Y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
        }else{
          qsad <- get(paste("q", sad, sep=""), mode = "function")
          ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(coef), dots))
        }
      }
    }else if(!missing(S) && !missing(A)){
      if (sad == "mand" || sad == "volkov" || sad== "mzsm" || sad == "ls" || sad == "geom" || sad == "nbinom"|| sad == "zipf"|| sad == "power"|| sad == "poilog"){
        y <- 1:A
        if(!missing(trunc)){
          X <- do.call(ptrunc, c(list(sad, q = y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
        }else{
          psad <- get(paste("p", sad, sep=""), mode = "function")
          X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(coef), dots))
        }
        f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
        ab <- f1(ppoints(S))
        #if(!any(is.na(ab[-1]))){
         # ab[1] <- sum(x) - sum(ab[-1])
        #}
      }else if(sad == "pareto" || sad == "gamma" || sad == "lnorm" || sad == "weibull"){
        Y <- ppoints(S)
        if(!missing(trunc)){
          ab <- do.call(qtrunc, c(list(sad, p = Y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
        }else{
          qsad <- get(paste("q", sad, sep=""), mode = "function")
          ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(coef), dots))
        }
      }
    }else 
      stop("inform S and A, or x")
  }else{
    S <- length(object@data$x)
    if (object@distr == "D"){
      y <- 1:max(object@data$x)
      if(!is.na(object@trunc)){
        if(object@sad=="ls") 
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(alpha = object@coef, N = sum(object@data$x)), lower.tail=F, trunc = object@trunc)))
        else if(object@sad=="zipf") 
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(s = object@coef, N = length(object@data$x)), lower.tail=F, trunc = object@trunc)))
        else if(object@sad=="mzsm") 
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(J = sum(object@data$x), theta = object@coef), lower.tail=F, trunc = object@trunc)))
        else if(object@sad=="volkov") 
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(J = sum(object@data$x), object@coef), lower.tail=F, trunc = object@trunc)))
        else if(object@sad=="mand") 
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(N = sum(object@data$x), object@coef), lower.tail=F, trunc = object@trunc)))
        else
          X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc)))
      }else {
        psad <- get(paste("p", object@sad, sep=""), mode = "function")
        if(object@sad=="ls")
          X <- do.call(psad, c(list(q = y, lower.tail = F, alpha = object@coef, N = sum(object@data$x))))
        else if(object@sad=="zipf")
          X <- do.call(psad, c(list(q = y, lower.tail = F, s = object@coef, N = length(object@data$x))))
        else if(object@sad=="mzsm")
          X <- do.call(psad, c(list(q = y, lower.tail = F, J = sum(object@data$x), theta = object@coef)))
        else if(object@sad=="volkov")
          X <- do.call(psad, c(list(q = y, lower.tail = F, J = sum(object@data$x)), as.list(object@coef)))
        else if(object@sad=="mand")
          X <- do.call(psad, c(list(q = y, lower.tail = F, N = sum(object@data$x)), as.list(object@coef)))
        else
          X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(object@coef)))
      }
      f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
      ab <- f1(ppoints(S))
      #if(is.na(ab[1]) & !any(is.na(ab[-1]))){
       # ab[1] <- sum(object@data$x) - sum(ab[-1])
      #}
    }else if(object@distr == "C"){
      Y <- ppoints(S)
      if(!is.na(object@trunc)){
        ab <- do.call(qtrunc, c(list(object@sad, p = Y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc)))
      }else{
        qsad <- get(paste("q", object@sad, sep=""), mode = "function")
        if(object@sad == "pareto")
          ab <- do.call(qsad, c(list(p = Y, lower.tail = F, scale = min(object@data$x)), as.list(object@coef)))
        else
          ab <- do.call(qsad, c(list(p = Y, lower.tail = F), as.list(object@coef)))
      }
    } else
      stop("unsupported distribution")
  }
  new("rad", data.frame(rank=1:S, abund=ab))
}
