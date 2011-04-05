dpoix <- function(y, frac, rate, log=FALSE) {  
	  b <- y*log(frac)
	  m <- log(rate)
	  n <- (y+1)*log(rate+frac)
    if(log)b+m-n else exp(b+m-n)
}
