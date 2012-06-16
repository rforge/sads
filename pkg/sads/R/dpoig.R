dpoig <- function(y, frac, rate, shape, trunc=0, log=FALSE) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
        abs(x - round(x)) < tol
    if(FALSE %in% sapply(y,is.wholenumber))
        stop("y must be integer because dpoig is a discrete PDF.")
    else {
        f <- function(y) {
            b <- y*log(frac)+shape*log(rate)+lgamma(y+shape)
            c <- lfactorial(y)+lgamma(shape)+(y+shape)*log(frac+rate)
            exp(b-c)
          }
        if(missing(trunc)) vals <- f(y)
        else vals <- f(y)/(1-f(trunc))
        if(log)log(vals) else vals
      }
}
