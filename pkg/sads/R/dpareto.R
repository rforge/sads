dpareto <- function(x, shape, scale = min(x), log = FALSE){
  if(shape <= 0 || scale <= 0)
    stop("shape and scale must be greater than zero")
  if(any(x < scale))
    stop("x must be greater than the scale parameter")
  lny <- log(shape) + shape*log(scale) - (shape+1)*log(x)
  if (log) return(lny)
  else return(exp(lny))
}