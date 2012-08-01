pls <- function(q, N, alpha, lower.tail=TRUE, log.p=FALSE){
  y <- sapply(q, function(x) sum(dls(1:x, N, alpha)))
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}