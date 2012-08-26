pzipf <- function(q, N, s){
  if (s <= 0) stop("s must be greater than zero")
  if (N < 1) stop("N must be positive integer")
  y <- NULL
  for (i in 1:length(q)){
    y[i] <- log(sum(1/(1:q[i])^s)) - log(sum(1/(1:N)^s))
  }
  return(exp(y))
}