dpoith <- function(y, frac, m=1, M=50, trunc=0) {
    f <- function(y) {
        b <- Igamma(y,frac*m,lower=FALSE) - Igamma(y,frac*M,lower=FALSE)
        c <- (factorial(y))*log(M/m)
        b/c
    }
    f(y)/(1-f(trunc))
}

