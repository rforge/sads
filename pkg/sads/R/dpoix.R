dpoix <- function(y, frac=0.05, rate=1/1000, log=FALSE) {  
	  b <- y*log(frac)
	  m <- log(rate)
	  n <- (y+1)*log(rate+frac)
    if(log)b+m-n else exp(b+m-n)
}