pmand <- function (q, N, s, v, lower.tail = TRUE, log.p = FALSE){
  if (s <= 0 ) 
    stop("s must be greater than zero")
  if (v <= 0) 
    stop("v must be greater than zero")
  if (N < 1||!any(is.wholenumber(N))) 
    stop("N must be positive integer")
  y <- NULL
  for (i in 1:length(q)) {
    y[i] <- log(sum(1/((1:q[i])+v)^s)) - log(sum(1/((1:N)+v)^s))
  }
  y <- exp(y)
  if (!lower.tail) 
    y <- 1 - y
  if (log.p) 
    y <- log(y)
  return(y)
}
