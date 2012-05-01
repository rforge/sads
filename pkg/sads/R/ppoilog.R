ppoilog <- function(q, mu, sig, lower.tail=TRUE, log.p=FALSE, trunc=0){
  if(!is.null(trunc)){
    if(any(q<=trunc))stop("at least one q larger than truncation value")
    f1 <- function(n){
      sum(dpoilog((trunc+1):n,mu,sig, trunc=trunc))
      }
    y <- sapply(q,f1)
  }
  else{
    y <- sapply(q,function(n)sum(dpoilog(0:n,mu,sig,trunc=trunc)))
  }
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}
