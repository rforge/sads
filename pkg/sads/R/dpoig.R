dpoig <- function(y, frac, rate, shape) {
  b <- y*log(frac)+shape*log(rate)+lgamma(y+shape)
  c <- lfactorial(y)+lgamma(shape)+(y+shape)*log(frac+rate)
  exp(b-c)
}
