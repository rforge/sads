integr2 <- function(n, J, m, theta, precisao = 0.01){
  res <- 0
  fun <- .C("Intgrl1", as.integer(n), as.integer(J), as.double(m),
            as.double(theta), as.double(precisao), result = res, PACKAGE = "zsm")
  fun[["result"]]
}

zsm <- function(n, J, m, theta, ...){
  f1 <- function(n){
    integr2(n, J, m, theta,...)
  }
  theta*sapply(n, f1)
}

dzsm <- function(n, J, m, theta, log = FALSE,...){
  all.values <- zsm(1:J, J, m, theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}