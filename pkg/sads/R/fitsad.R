fitsad <- function(x, sad =c("gamma","geom","lnorm","ls","mzsm","nbinom","pareto","poilog","power","weibull","volkov"), ...){ 
  dots <- list(...)
  fit <- get(paste("fit", sad, sep=""), mode = "function")
  if(!"trunc" %in% names(dots) & sad=="poilog")
    do.call(fit, c(list(x = x, trunc = 0), dots))
  else
    do.call(fit, c(list(x = x), dots))
}
