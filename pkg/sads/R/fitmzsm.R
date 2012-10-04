fitmzsm <- function(x, trunc, start.value, ...){
  dots <- list(...)
  J <- sum(x)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    theta <- mean(x)
  } else{
    theta <- start.value
  }
  if (missing(trunc)){
    LL <- function(theta) -sum(dmzsm(x, J, theta, log = TRUE))
  } else {
    LL <- function(theta) -sum(dtrunc("mzsm", x, coef = list(J = J, theta = theta), trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(theta = theta), data = list(x = x), ...)
  new("fitsad", result, sad="mzsm", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc)) 
}