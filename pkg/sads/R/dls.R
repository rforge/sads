dls <- function(x, N, alpha, log = FALSE){
  if (!all(is.finite(c(N, alpha)))) stop("all parameters should be finite")
  if (N <= 0)  stop("N must be larger than zero")
  if (alpha <= 0)  stop("alpha must be larger than zero")
  if(any(x < 1)) stop("at least one x less than one")
  if(any(!is.wholenumber(x))) stop("at least one non-integer x or all x must be integers")
  X <- N/(N+alpha)
  gama <- function(y) {1/log(1/(1-y))}
  y <- gama(X)*(X^x)/x
  if(log) return (log(y))
  else return(y)
}