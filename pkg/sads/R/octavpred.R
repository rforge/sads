octavpred <- function(object, x, sad, coef, trunc, oct, S,...){
  dots <- list(...)
  if (!missing(sad) && !missing(coef)){
    if(!missing(x)){
      S <- length(x)
      if(missing(oct)){
        oct <- 1:(ceiling(max(log2(x)))+1)
        if(any(x < 1)){
          octlower <- ceiling(min(log2((x)))+1):0
          oct <- c(octlower, oct)
        }
      }
    }else if(!missing(oct)){
      oct <- 1:oct
    }else stop("inform x or oct")
    n <- 2^(oct-1)
    if(!missing(trunc)){
        Y <- do.call(ptrunc, c(list(sad, q = n, coef = as.list(coef), trunc = trunc), dots))
    }else {
      psad <- get(paste("p", sad, sep=""), mode = "function")
      Y <- do.call(psad, c(list(q = n), as.list(coef), dots))
    }
  }else{
    S <- length(object@data$x)
    if(missing(oct)){
      oct <- 1:(ceiling(max(log2(object@data$x)))+1)
      if(any(object@data$x < 1)){
        octlower <- ceiling(min(log2((object@data$x)))+1):0
        oct <- c(octlower, oct)
      }
    }
    n <- 2^(oct-1)
    if(!is.na(object@trunc)){
      if(object@sad == "ls")
        Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = list(alpha = object@coef, N = sum(object@data$x)), trunc = object@trunc), dots))
      else
        Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = as.list(object@coef), trunc = object@trunc), dots))
    }else {
      psad <- get(paste("p", object@sad, sep=""), mode = "function")
      if(object@sad == "ls")
        Y <- do.call(psad, c(list(q = n, alpha = object@coef, N = sum(object@data$x)), dots))
      else
        Y <- do.call(psad, c(list(q = n), as.list(object@coef), dots))
    }
  }
  Y <- c(Y[1], diff(Y))*S
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}
