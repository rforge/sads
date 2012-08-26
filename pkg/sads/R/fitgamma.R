fitgamma <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    ka <- mean(x)/sd(x)
    theta <- ka^2
  } else{
    ka <- start.value[1]
    theta <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(shape, scale) -sum(dgamma(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(trunc("dgamma", x, shape, scale, trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(shape = ka, scale = theta), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="gamma", trunc = ifelse(missing(trunc), NaN, trunc)) 
}