require("zipfR")
dpoith <- function(y, frac, m=0.1, M=50000) {
  b <- Igamma(y,frac*m,lower=FALSE) - Igamma(y,frac*M,lower=FALSE)
  c <- factorial(y)*log(M/m)
  b/c
}
