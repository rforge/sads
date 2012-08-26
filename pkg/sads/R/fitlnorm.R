fitlnorm <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    meanlog <- sum(log(x))/length(x)
    sdlog <- sqrt((sum(log(x)-meanlog)^2)/length(x))
  } else{
    meanlog <- start.value[1]
    sdlog <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(meanlog, sdlog) -sum(dlnorm(x, meanlog, sdlog, log = TRUE))
  } else {
    LL <- function(meanlog, sdlog) -sum(trunc("dlnorm", x, meanlog, sdlog, trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(meanlog = meanlog, sdlog = sdlog), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="lnorm", trunc = ifelse(missing(trunc), NaN, trunc)) 
}