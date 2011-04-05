library("zipfR")
dpoith <- function(y, frac, m=100, M=5000) {
  b <- Igamma(y,frac*m,lower=FALSE) - Igamma(y,frac*M,lower=FALSE)
  c <- factorial(y)*log(M/m)
  b/c
}
