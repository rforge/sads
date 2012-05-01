octavpred <- function(x,rich,sad,oct,coef,...){
  if(missing(rich)) S <- length(x)
  if(missing(oct)&!missing(x))oct <- 1:(ceiling(max(log2((x))))+1)
  n <- 2^(oct-1)
  psad <- paste("p",sad,sep="")
  if(missing(coef))dots <- list(...)
  else dots <- coef
  if(sad=="ls"&!"N"%in%names(dots)) dots$N <- sum(x)
  Y <- do.call(psad,c(list(q=n),dots))
  Y <- c(Y[1],diff(Y))*S
  new("octav",data.frame(octave=oct, upper=factor(n), Freq=Y))
}
