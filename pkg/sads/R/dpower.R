dpower <- function(x, s, log = FALSE){
  if (any(x < 1)) warning("the zipf's distribution is not set to zero")
  if (s <= 1) stop("s must be greater than one")
  if (!any(is.wholenumber(x))) warning("x must be integer")
  y <- NULL
  for (i in 1:length(x)){
    if(!is.wholenumber(x[i])) y[i] <- -Inf
    else y[i] <- -s*log(x[i])-log(zeta(s))
  }
  if(log) return(y)
  else return(exp(y))
}
