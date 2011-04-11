library(fAsianOptions)
dnex <- function(y, frac, rate, k=0.5) {
    b <- rate*k*kummerM(1+y,2-k,k*rate/frac)
    c <- frac*(k-1)
	  m <- b/c
    
    b <- ((k*rate)^k)*pi*(1/sin(k*pi))*gamma(k+y)*kummerM(k+y,k,k*rate/frac)
    c <- (frac^k)*(gamma(k)^2)*gamma(y+1)
	  n <- b/c
    
    m+n
}
