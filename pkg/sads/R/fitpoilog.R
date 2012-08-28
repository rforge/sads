fitpoilog <- function(x, trunc, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(trunc)){
    pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))), zTrunc = F)$par
    LL <- function(mu, sig) -sum(dpoilog(x, mu, sig, log = TRUE))
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), ...)
  } else{
    pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
    LL <- function(mu, sig) -sum(trunc("dpoilog", x, mu, sig, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), ...)
  }
  new("fitsad", result, sad="poilog", trunc = ifelse(missing(trunc), NaN, trunc))
}