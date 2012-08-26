fitnbinom <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    LL <- function(prob) -sum(dnbinom(x, size = length(x), prob, log = TRUE))
    phat<- length(x)/(length(x)+mean(x))
    phat <- mle2(LL, start = list(prob = phat), method="Brent", upper = 1, lower = 0)@coef
    mu <- length(x)*(1 - phat)/phat
  } else{
    mu <- start.value
  }
  if (missing(trunc)){
    LL <- function(mu) -sum(dnbinom(x, size = length(x), mu, log = TRUE))
  } else{
    LL <- function(mu) -sum(trunc("dnbinom", x, size = length(x), mu, trunc = trunc, log = TRUE))
  }
  result <- mle2(LL, start = list(mu = mu), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="nbinom", trunc = ifelse(missing(trunc), NaN, trunc))
}