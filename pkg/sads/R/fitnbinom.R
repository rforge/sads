fitnbinom <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    phat <- length(x)/(length(x) + mean(x))
  } else{
    phat <- start.value
  }
  if (missing(trunc)){
    LL <- function(prob) -sum(dnbinom(x, size = length(x), prob, log = TRUE))
  } else{
    LL <- function(prob) -sum(dtrunc("nbinom", x = x, coef = list(size = length(x), prob), trunc = trunc, log = TRUE))
  }
  result <- mle2(LL, start = list(prob = phat), data = list(x = x), method="Brent", upper = 1, lower = 0, ...)
  new("fitsad", result, sad="nbinom", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc))
}