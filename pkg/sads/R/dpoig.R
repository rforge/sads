dpoig <- function(y, frac, rate, shape, trunc=0) {
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
        abs(x - round(x)) < tol
    if(FALSE %in% sapply(y,is.wholenumber))
        print("y must be integer because dpoig is a discrete PDF.")
    else {
        f <- function(y) {
            b <- y*log(frac)+shape*log(rate)+lgamma(y+shape)
            c <- lfactorial(y)+lgamma(shape)+(y+shape)*log(frac+rate)
            exp(b-c)
        }
        f(y)/(1-f(trunc))
    }
}
