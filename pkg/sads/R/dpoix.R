dpoix <- function(x, frac, rate, log=FALSE) {
  
	  b <- x*log(frac)
	  m <- log(rate)
	  n <- (x+1)*log(rate+frac)
    if(log)b+m-n else exp(b+m-n)
}
