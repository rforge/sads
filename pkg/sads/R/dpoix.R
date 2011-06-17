dpoix <- function(y, frac, rate, log=FALSE) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
        abs(x - round(x)) < tol
    if(FALSE %in% sapply(y,is.wholenumber))
       print("y must be integer because dpoix is a discrete PDF.")
    else {
        b <- y*log(frac)
	      m <- log(rate)
	      n <- (y+1)*log(rate+frac)
        if(log)b+m-n else exp(b+m-n)
    }
}