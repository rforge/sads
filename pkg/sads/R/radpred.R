radpred <- function(x, sad, rad, coef, trunc, S, A, N, ...){
  dots <- list(...)
  if(missing(trunc)) trunc <- NA
  if(class(x)=="fitrad"||!missing(rad)){
    if(!missing(x)){
      if (class(x)=="numeric"){
        S <- length(x)
        N <- sum(x)
      }
      else if (class(x)=="fitrad"){
        S <- length(x@rad.tab$abund)
        N <- sum(x@rad.tab$abund)
        coef <- bbmle::coef(x)
        rad <- x@rad
        trunc <- x@trunc
      }
    }
    else if(missing(S) || missing(N)) stop("please provide x or S and N")
        
    if (rad == "zipf" || rad == "mand"){
      y <- 1:S
      if(!is.na(trunc)){
        ab <- do.call(dtrunc, c(list(rad, x = y, coef = as.list(coef), trunc = trunc), dots))*N
      }
      else{
        drad <- get(paste("d", rad, sep=""),  mode = "function")
        ab <- do.call(drad, c(list(x = y), as.list(coef), dots))*N
      }
    }
    else stop("unsupported distribution")
  } 
      
  else if(class(x)=="fitsad"||!missing(sad)){
    if(!missing(sad) && !missing(coef)){
      if (!missing(x) && class(x)=="numeric"){
        S <- length(x)
        if (sad == "volkov" || sad == "mzsm" || sad == "ls" || sad == "geom" || sad == "nbinom"|| sad == "power"|| sad == "poilog"){
          y <- 1:max(x)
          if(!is.na(trunc)){
            X <- do.call(ptrunc, c(list(sad, q = y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
          }
          else{
            psad <- get(paste("p", sad, sep=""),  mode = "function")
            X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(coef), dots))
          }
          f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
          ab <- f1(ppoints(S))
          ##if(!any(is.na(ab[-1]))){
          ##  ab[1] <- sum(x) - sum(ab[-1])
          ##}
        }
        else if(sad == "pareto" || sad == "gamma" || sad == "lnorm" || sad == "weibull"){
          Y <- ppoints(S)
          if(!is.na(trunc)){
            ab <- do.call(qtrunc, c(list(sad, p = Y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
          }
          else{
            qsad <- get(paste("q", sad, sep=""), mode = "function")
            ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(coef), dots))
          }
        }
      }
      else if(!missing(S) && !missing(A)){
        if (sad == "volkov" || sad== "mzsm" || sad == "ls" || sad == "geom" || sad == "nbinom"||sad == "power"|| sad == "poilog"){
          y <- 1:A
          if(!is.na(trunc)){
            X <- do.call(ptrunc, c(list(sad, q = y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
          }
          else{
            psad <- get(paste("p", sad, sep=""), mode = "function")
            X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(coef), dots))
          }
          f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
          ab <- f1(ppoints(S))
          ##if(!any(is.na(ab[-1]))){
          ## ab[1] <- sum(x) - sum(ab[-1])
          ##}
        }
        else if(sad == "pareto" || sad == "gamma" || sad == "lnorm" || sad == "weibull"){
          Y <- ppoints(S)
          if(!is.na(trunc)){
            ab <- do.call(qtrunc, c(list(sad, p = Y, coef = as.list(coef), lower.tail=F, trunc = trunc), dots))
          }
          else{
            qsad <- get(paste("q", sad, sep=""), mode = "function")
            ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(coef), dots))
          }
        }
      }
      else 
        stop("input S and A, or x")
    }
    else{
      S <- length(x@data$x)
      if (x@distr == "D"){
        y <- 1:max(x@data$x)
        if(!is.na(x@trunc)){
          if(x@sad=="ls") 
            X <- do.call(ptrunc, c(list(x@sad, q = y, coef = list(alpha = x@coef, N = sum(x@data$x)),
                                        lower.tail=F, trunc = x@trunc)))
          else if(x@sad=="mzsm") 
            X <- do.call(ptrunc, c(list(x@sad, q = y, coef = list(J = sum(x@data$x), theta = x@coef),
                                        lower.tail=F, trunc = x@trunc)))
          else if(x@sad=="volkov") 
            X <- do.call(ptrunc, c(list(x@sad, q = y, coef = as.list(J = sum(x@data$x), x@coef),
                                        lower.tail=F, trunc = x@trunc)))
          else
            X <- do.call(ptrunc, c(list(x@sad, q = y, coef = as.list(x@coef), lower.tail=F, trunc = x@trunc)))
        }
        else {
          psad <- get(paste("p", x@sad, sep=""), mode = "function")
          if(x@sad=="ls")
            X <- do.call(psad, c(list(q = y, lower.tail = F, alpha = x@coef, N = sum(x@data$x))))
          else if(x@sad=="mzsm")
            X <- do.call(psad, c(list(q = y, lower.tail = F, J = sum(x@data$x), theta = x@coef)))
          else if(x@sad=="volkov")
            X <- do.call(psad, c(list(q = y, lower.tail = F, J = sum(x@data$x)), as.list(x@coef)))
          else
            X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(x@coef)))
        }
        f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
        ab <- f1(ppoints(S))
        ##if(is.na(ab[1]) & !any(is.na(ab[-1]))){
        ## ab[1] <- sum(x@data$x) - sum(ab[-1])
        ##}
      }
      else if(x@distr == "C"){
        Y <- ppoints(S)
        if(!is.na(x@trunc)){
          ab <- do.call(qtrunc, c(list(x@sad, p = Y, coef = as.list(x@coef), lower.tail=F, trunc = x@trunc)))
        }
        else{
          qsad <- get(paste("q", x@sad, sep=""), mode = "function")
          if(x@sad == "pareto")
            ab <- do.call(qsad, c(list(p = Y, lower.tail = F, scale = min(x@data$x)), as.list(x@coef)))
          else
            ab <- do.call(qsad, c(list(p = Y, lower.tail = F), as.list(x@coef)))
        }
      }
      else
        stop("unsupported distribution")
    }
  }
  else if(missing(sad) && missing(rad)) stop ("please input sad OR rad") 
  new("rad", data.frame(rank=1:S, abund=ab))
}
