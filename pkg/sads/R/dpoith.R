require("zipfR")
dpoith <- function(y, frac, m=1, M=50) {
  b <- Igamma(y,frac*m,lower=FALSE) - Igamma(y,frac*M,lower=FALSE)
  c <- (factorial(y))*log(M/m)
  b/c
}