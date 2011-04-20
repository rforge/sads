dpoig <- function(y, frac, rate, shape) {
  b <- (frac^y)*(rate^shape)*gamma(y+shape)
  c <- factorial(y)*gamma(shape)*(frac+rate)^(y+shape)
  b/c
}