dls <- function(x, N, alpha, log = FALSE){
  if(any(x < 1)) stop("at least one x less than one")
  if(any(!is.wholenumber(x))) stop("at least one non-integer x")
  X <- N/(N+alpha)
  gama <- function(y) {1/log(1/(1-y))}
  y <- gama(X)*(X^x)/x
  if(log) log(y)
  else y
}
