dpoix <- function(y, frac, rate, trunc=0, log=FALSE) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
      abs(x - round(x)) < tol
    }
    if(sum(y,is.wholenumber(y))<length(y))
       stop("dpoix is a discrete PDF; all y's must be integers")
    else {
      f <- function(y){
        b <- y*log(frac)
        m <- log(rate)
        n <- (y+1)*log(rate+frac)
        exp(b+m-n)
      }
      if(!is.null(trunc)) samp <- f(y)/(1-sum(f(0:trunc))) else samp <- f(y)
      if(log)samp else log(samp)
    }
  }
       
    
