ppoilog2 <- function(q, mu, sig, lower.tail=TRUE, log.p=FALSE){
  y <- sapply(q, function(n) sum(poilog::dpoilog(0:n, mu, sig)))
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}