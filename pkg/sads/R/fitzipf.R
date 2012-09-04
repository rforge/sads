fitzipf <- function(x, N, trunc, start.value, upper = 20, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(N)){
    N <- length(x)
  }
  if(missing(start.value)){
    p <- x/sum(x)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    sss <- opt$minimum
  }else{
    sss <- start.value
  }
  if(missing(trunc)){
    LL <- function(s) -sum(dzipf(x, N, s, log = TRUE))
  } else{
    LL <- function(s) -sum(dtrunc("zipf", x = x, coef = list(N = N, s = s), trunc = trunc, log = TRUE))
  }
  result <-  mle2(LL, start = list(s = sss), data = list(x = x), method = "Brent", lower = 0, upper = upper, ...)
  if(abs(as.numeric(result@coef) - upper) < 0.001) warning("mle equal to upper bound provided. \n Try value for the 'upper' arguent")
  new("fitsad", result, sad="zipf", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc))
}
