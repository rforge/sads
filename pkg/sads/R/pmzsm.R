pmzsm <- function(q, J, theta, lower.tail=TRUE, log.p=FALSE){
  if (length(J) > 1 | length(theta) > 1) stop("vectorization of J and theta is currently not implemented")
  if (!all(is.finite(c(J, theta)))) stop("all parameters should be finite")
  if (J <= 0)  stop("J must be larger than zero")
  if (theta <= 0) stop("theta must be larger than zero")
  z <- cumsum(dmzsm(1:max(q), J, theta))
  y <- z[q]
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}