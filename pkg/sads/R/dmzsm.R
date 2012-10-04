mzsm <- function(n, J, theta){
  lengn <- length(n)
  res <- rep(0, lengn)
  fun <- .C("msn", as.integer(n), as.integer(lengn), as.integer(J), 
            as.double(theta), result = res, PACKAGE = "zsm")
  fun[["result"]]
}

dmzsm <- function(n, J, theta, log = FALSE,...){
  all.values <- mzsm(1:J, J, theta)
  all.values[all.values == Inf] <- NaN
  lprobs <- log(all.values[n])-log(sum(all.values, na.rm=TRUE))
  if(log) lprobs
  else exp(lprobs)
}