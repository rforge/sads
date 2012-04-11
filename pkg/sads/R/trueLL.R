trueLL <- function(x,dens,precision=1,log=TRUE,...){
  cdf <- paste("p",deparse(substitute(dens)),sep="")
  dots <- list(...)
  D <- precision/2
  probs <- do.call(cdf, c(list(q=x+D),dots))-
    do.call(cdf, c(list(q=x-D),dots))
  if(log)sum(log(probs)) else prod(probs)
}
