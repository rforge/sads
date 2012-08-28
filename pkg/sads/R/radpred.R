radpred <- function(x,sad,coef,...){
  psad <- paste("p",sad,sep="")
  if(missing(coef))dots <- list(...)
  else dots <- c(coef, list(...))
  if(sad=="ls"&!"N"%in%names(dots)) dots$N <- sum(x)
  S <- length(x)
  y <- 1:max(x)
  X <- do.call(psad,c(list(q=y,lower.tail=F),dots))
  f1 <- approxfun(x=c(1,X),y=c(0,y), method="constant")
  ab <- f1(ppoints(S))
  if(is.na(ab[1])&!any(is.na(ab[-1]))){
    ab[1] <- sum(x)-sum(ab[-1])
  }
  new("rad",data.frame(rank=1:S, abund=ab))
}
