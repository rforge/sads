fitpoilog2 <- function(x, trunc = 0, ...){
  if (trunc == 0){
    result <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))
  } else{
    pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
    LL <- function(mu, sig) -sum(trunc("dpoilog", x, mu, sig, trunc, log = TRUE))
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), method = "SANN", ...)
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), ...)
  }
  new("fitsad", result, sad="poilog", trunc = trunc)
}
