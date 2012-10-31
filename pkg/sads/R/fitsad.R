fitsad <- function(x, sad, trunc, ...){ #trueLL dec.places - continuas
  dots <- list(...)
  fit <- get(paste("fit", sad, sep=""), mode = "function")
  if(sad=="volkov" || sad =="mzsm" || sad == "geom" || sad == "ls" || sad == "nbinom" || sad == "poilog" || sad == "power" || sad == "zipf" || sad == "mand"){
    if(missing(trunc)){
      if(sad == "poilog") {
        trunc = 0
        do.call(fit, c(list(x = x), trunc = trunc, dots))
      }else
        do.call(fit, c(list(x = x), dots))
    }else {
      do.call(fit, c(list(x = x), trunc = trunc, dots))
    }
  }else if(sad == "pareto" || sad == "gamma" || sad == "lnorm" || sad == "weibull"){
    if(missing(trunc)){
      do.call(fit, c(list(x = x), dots))
    }else{
      do.call(fit, c(list(x = x), trunc = trunc, dots))
    }
  }else{
    stop("unsupported distribution")
  }
}
