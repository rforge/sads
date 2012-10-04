fitsad <- function(x, sad, trunc, ...){ #trueLL dec.places - continuas
  dots <- list(...)
  fit <- get(paste("fit", sad, sep=""), mode = "function")
  if(sad=="zsm" || sad =="mzsm" || sad == "geom" || sad == "ls" || sad == "nbinom" || sad == "poilog" || sad == "power" || sad == "zipf"){
    if(missing(trunc)){
      if(sad == "geom" || sad == "nbinom" || sad == "poilog") {
        trunc = 0
        do.call(fit, c(list(x = x), trunc = trunc, dots))
      }else
        do.call(fit, c(list(x = x), dots))
    }else {
      do.call(fit, c(list(x = x), trunc = trunc, dots))
    }
  }else if(sad == "gamma" || sad == "lnorm" || sad == "weibull"){
    if(missing(trunc)){
      do.call(fit, c(list(x = x), dots))
    }else{
      do.call(fit, c(list(x = x), trunc = trunc, dots))
    }
  }else{
    stop("unsupported distribution")
  }
}
