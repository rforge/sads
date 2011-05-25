dneg <- function(y, frac, rate, k=2, beta=1) {
  s <- exp(lgamma(k+y)-(lgamma(k)+lgamma(y)))
  b <- pi*s*(1/sin(pi*(beta-k)))
  c <- gamma(beta)*gamma(k+y)

  t1 <- (k*rate/frac)^k
  t2 <- gamma(k+y)
  t3 <- kummerM(k+y,1-beta+k,k*rate/frac)
  d1 <- t1*t2*t3
  
  t4 <- gamma(beta+y)
  t5 <- kummerM(beta+y,1+beta-k,k*rate/frac)
  d2 <- t1*t4*t5
  
  b/c*(d1-d2)
}
