fitzipf <- function(x, N, trunc = 0, start.value, upper = 20, ...){
  if (min(x)<=trunc){
    stop("truncation point should be lower than the lowest data value")
  }
  if(missing(N)){
    N <- length(x)
  }
  if(missing(start.value)){
    p <- x/sum(x)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    print(opt)
    sss <- opt$minimum
    #LL1 <-function(s) sum(x*(s*log(1:N)+log(sum(1/(1:N)^s))))
    #fit <- mle2(LL1, start = list(s = sss))
    #sss <- exp(fit@coef)
    #print(fit)
  }else{
    sss <- start.value
  }
  if(trunc >= 1){
    LL <- function(s) -sum(trunc("dzipf", x, N, s, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = list(s = sss), data = list(x = x), method = "Brent", lower = 0, upper = upper, ...)
  } else{
    LL <-function(s) sum(x*(((1:N)^s)*(sum(1/(1:N)^s))))
    result <- mle2(LL, start = list(s = sss), data = list(x = x), method="Brent", lower = 0, upper = upper, ...)
    #LL <- function(s) -sum(dzipf(x, N, s, log = TRUE))
    #result <- mle2(LL, start = list(s = sss), data = list(x = x), method="Brent", lower = 0, upper = upper, ...)
  }
  if(abs(as.numeric(result@coef) - upper) < 0.001) warning("Check the upper limit of the function")
  new("fitsad", result, sad="zipf", trunc = trunc)
}