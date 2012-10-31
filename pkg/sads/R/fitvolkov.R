fitvolkov <- function(x, trunc, start.value,...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    thetahat <- mean(x)
    mhat <- runif(1)
  } else{
    mhat <-start.value[1]
    thetahat <- start.value[2]
  }
  if (missing(trunc)){
    LL <- function(m, theta) -sum(dvolkov(x, J = sum(x), m = m, theta = theta, log = TRUE))
  } else {
    LL <- function(m, theta) -sum(dtrunc("volkov", x = x, coef = list(c(J = sum(x), m = m, theta = theta)), trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(m = mhat, theta = thetahat), data = list(x = x), method="L-BFGS-B", lower=c(1e-4, 1), upper=c(1, length(x)))
  new("fitsad", result, sad="volkov", distr = "D", trunc = ifelse(missing(trunc), NaN, trunc)) 
}