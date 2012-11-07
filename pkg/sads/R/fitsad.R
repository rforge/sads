fitsad <- function(x, sad =c("gamma","geom","lnorm","ls","mzsm","nbinom","pareto","poilog","power","weibull","volkov"),
trunc, start.value, trueLL = TRUE, dec.places = 0, ...){ 
  dots <- list(...)
  fit <- get(paste("fit", sad, sep=""), mode = "function")
  if(sad=="volkov" || sad =="mzsm" || sad == "geom" || sad == "ls" || sad == "nbinom" || sad == "poilog" || sad == "power"){
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
