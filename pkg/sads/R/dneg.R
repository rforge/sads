dneg <- function(y, frac, rate, shape, k) {
  s <- exp(lgamma(k+y)-(lgamma(k)+lgamma(y)))
  b <- pi*s*(1/sin(pi*(shape-k)))
  c <- gamma(shape)*gamma(k+y)

  t1 <- (k*rate/frac)^k
  t2 <- gamma(k+y)
  t3 <- kummerM(k+y,1-shape+k,k*rate/frac)
  d1 <- t1*t2*t3
  
  t4 <- gamma(shape+y)
  t5 <- kummerM(shape+y,1+shape-k,k*rate/frac)
  d2 <- t1*t4*t5
  
  b/c*(d1-d2)
}
