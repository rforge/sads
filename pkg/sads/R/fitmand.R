fitmand <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    shat <- 2
    vhat <- 30
  } else{
    shat <- start.value[1]
    vhat <- start.value[2]
  }
  if (missing(trunc)){
    LL <- function(s, v) -sum(dmand(x, N = sum(x), s, v, log = TRUE))
  } else{
    LL <- function(s, v) -sum(dtrunc("mand", x = x, coef = list(N = sum(x), s = s, v = v), trunc = trunc, log = TRUE))
  }
  result <- mle2(LL, start = list(s = shat, v = vhat), data = list(x = x), ...)
  new("fitsad", result, sad="mand", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc))
}