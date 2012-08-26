ppower <- function(q, s){
  if (s <= 1) stop("s must be greater than one")
  y <- NULL
  for (i in 1:length(q)){
    y[i] <- log(sum(1/(1:q[i])^s)) - log(zeta(s))
  }
  return(exp(y))
}