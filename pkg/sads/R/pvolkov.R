pvolkov <- function(q, theta, m , J, lower.tail=TRUE, log.p=FALSE, ...){
  if(length(J) > 1 | length(theta) > 1) stop("vectorization of J and theta is currently not implemented")
  if(length(m)>1) stop("vectorization of m is currently not implemented")
  if(!all(is.finite(c(J, theta, m)))) stop("all parameters should be finite")
  if(J <= 0)  stop("J must be larger than zero")
  if(theta <= 0) stop("theta must be larger than zero")
  if(m < 0 | m > 1) stop("m must be between zero and one")
  z <- cumsum(dvolkov(x = 1:max(q), theta = theta, J = J, m = m, ...))
  y <- z[q]
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}
