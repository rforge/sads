dpoig <- function(y, prob, rate=1, size) {
  frac <- (1 - prob)/prob
  shape <- size
  b <- (frac^y)*(rate^shape)*gamma(y+shape)
  c <- factorial(y)*gamma(shape)*(frac+rate)^(y+shape)
  b/c
}
