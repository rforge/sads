trueLL <- function(x,dens,precision=1,trunc,log=TRUE,...){
  cdf <- paste("p",deparse(substitute(dens)),sep="") 
  dots <- list(...)
  D <- precision/2
  if(missing(trunc)){
    probs <- do.call(cdf, c(list(q=x+D),dots))-
    do.call(cdf, c(list(q=x-D),dots))
  }
  else{
    probs <- do.call("trunc",c(list(f=cdf, x=x+D, trunc=trunc),dots))-
    do.call("trunc", c(list(f=cdf,q=x-D, trunc=trunc),dots))
  }
  if(log)sum(log(probs)) else prod(probs)
}
