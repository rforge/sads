library("zipfR")
dpoith <- function(x, frac, m=100, M=5000) {
  b <- Igamma(x,frac*m,lower=FALSE) - Igamma(x,frac*M,lower=FALSE)
  c <- factorial(x)*log(M/m)
  b/c
}
