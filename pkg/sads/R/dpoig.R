dpoig <- function(x, frac, rate, shape) {
  
  b <- (frac^x)*(rate^shape)*gamma(x+shape)
  c <- factorial(x)*gamma(shape)*(frac+rate)^(x+shape)
  b/c
}
