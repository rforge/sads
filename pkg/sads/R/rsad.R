rsad <- function(S,frac,sad, Pois.samp=TRUE,k,zeroes=FALSE,...){
  sad <- paste("r",sad,sep="")
  dots <- list(...)
  com <- do.call(sad,c(list(n=S),dots))
  if(Pois.samp) sam=rpois(S,lambda=frac*com)
  else {
      if(missing(k))stop("For negative binomial sampling please provide a value for k")
      else sam=rnbinom(S,mu=frac*com,size=k)}
  if(zeroes)sam else sam[sam>0]
}
