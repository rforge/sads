dneth <- function(y, frac, rate, k=2, M=1, m=10) {
   d <- 1/(y*log(M/m))
   b <- (k^k)*gamma(k+y)
   c <- gamma(y+1)*(k-1)
   
   F2 <- kummerU(k-1,k+y,k,-k/frac*m)
   d1 <- m*F2/(frac*m)^k
   
   F2 <- kummerU(k-1,k+y,k,-k/frac*M)
   d2 <- M*F2/(frac*M)^k
   
   d*(b/c*(d1-d2))
}
