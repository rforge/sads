fitzsm <- function(x, trunc, start.value,...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    theta <- mean(x)
    m <- runif(1)
  } else{
    m <-start.value[1]
    theta <- start.value[2]
  }
  J <- sum(x)
  if (missing(trunc)){
    LL <- function(m, theta) -sum(dzsm(x, J, m, theta, log = TRUE))
  } else {
    LL <- function(m, theta) -sum(dtrunc("zsm", x, coef = list(c(J = J, m = m, theta = theta)), trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list( m = m, theta = theta), data = list(x = x), method="L-BFGS-B", lower=c(0.00001, 1), upper=c(0.9999999, length(x)))
  new("fitsad", result, sad="zsm", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc)) 
}