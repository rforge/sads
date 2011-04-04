dpoig <- function(x, prob, rate=1, size) {
  frac <- (1 - prob)/prob
  shape <- size
  b <- (frac^x)*(rate^shape)*gamma(x+shape)
  c <- factorial(x)*gamma(shape)*(frac+rate)^(x+shape)
  b/c
}
