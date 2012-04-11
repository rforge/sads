rsad <- function(S,frac,sad, samp=c("Poisson","Negbinom"),k,zeroes=FALSE,...){
  sad <- paste("r",deparse(substitute(sad)),sep="")
  dots <- list(...)
  com <- do.call(sad,c(list(n=S),dots))
  sam <- switch(samp,
                Poisson=rpois(S,lambda=frac*com),
                Negbinom=rnbinom(S,mu=frac*com,size=k)
                )
  if(zeroes)sam else sam[sam>0]
}
