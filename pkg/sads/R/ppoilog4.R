ppoilog4 <- function(q, mu, sig, lower.tail=TRUE, log.p=FALSE){
  z <- cumsum(poilog::dpoilog(0:max(q),mu,sig))
  y <- z[q+1]
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}
