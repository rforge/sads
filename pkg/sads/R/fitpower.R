fitpower <- function(x, trunc, start.value, upper = 20, ...){
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    shat <- 2
  } else{
    shat <- start.value
  }
  if (missing(trunc)){
    LL <- function(s) -sum(dpower(x, s, log = T))
  } else{
    LL <- function(s) -sum(trunc("dpower", x, s, trunc = trunc, log = T))
  }
  result <- mle2(LL, start = list(s = shat), data = list(x = x), method = "Brent", lower = 1, upper = upper, ...)
  if(abs(as.numeric(result@coef) - upper) < 0.0000001) warning("mle equal to upper bound provided. \n Try value for the 'upper' arguent")
  new("fitsad", result, sad = "power", trunc = ifelse(missing(trunc), NaN, trunc))
}
