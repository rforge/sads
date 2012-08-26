fitgamma <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    ka <- (mean(x)/sd(x))^2
    theta <- var(x)/mean(x)
    kahat <- function(k, dados){
      eq <- length(dados)*(log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
      eq
    }
    ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x)$root
    theta <- mean(x)/ka
  } else{
    ka <- start.value[1]
    theta <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(shape, scale) -sum(dgamma(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(trunc("dgamma", x, shape, scale, trunc = trunc, log = TRUE))
  }  
  #result <- mle2(LL, start = list(shape = ka, scale = theta), method="SANN")
  result <- mle2(LL, start = list(shape = ka, scale = theta), data = list(x = x), ...)
  new("fitsad", result, sad="gamma", trunc = ifelse(missing(trunc), NaN, trunc)) 
}