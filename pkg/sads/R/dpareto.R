dpareto <- function(x, shape, scale = min(x), log = FALSE){
  if(shape <= 0 || scale <= 0)
    stop("shape and scale must be greater than zero")
  if(any(x < scale))
    stop("scale parameter must be equal or greater than min(x)")
  lny <- log(shape) + shape*log(scale) - (shape+1)*log(x)
  if (log) return(lny)
  else return(exp(lny))
}
