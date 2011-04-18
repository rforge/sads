require(gtools)
require(fAsianOptions)
dneg <- function(y, frac, rate, k=2, beta=1) {
    s <- length(combinations(k+y-1,k-1)[,1])
    b <- pi*s*(1/sin(pi*(beta-k)))
    c <- gamma(beta)*gamma(k+y)
	  
    t1 <- (k*rate/frac)^k
    t2 <- gamma(k+y)
    t3 <- kummerM(k+y,1-beta+k,k*rate/frac)
    d1 <- t1*t2*t3
    
    t1 <- (k*rate/frac)^k
    t2 <- gamma(beta+y)
    t3 <- kummerM(beta+y,1+beta-k,k*rate/frac)
    d2 <- t1*t2*t3
    
    b/c*(d1-d2)
}
