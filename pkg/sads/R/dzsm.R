integr2 <- function(n, J, m, theta, precisao = 0.01){
  res <- 0
  fun <- .C("Intgrl1", as.integer(n), as.integer(J), as.double(m),
            as.double(theta), as.double(precisao), result = res, PACKAGE = "sads")
  fun[["result"]]
}

zsm <- function(n, J, m, theta, ...){
  f1 <- function(n){
    integr2(n, J, m, theta,...)
  }
  theta*sapply(n, f1)
}

dzsm <- function(x, J, m, theta, log = FALSE){
  if (J <= 0) stop ("J must be great than zero")
  if (theta <= 0) stop ("theta must be great than zero")
  if (m > 1 | m < 0) stop ("m must be between zero and one")
  sn <- zsm(y=x, J = J, m = m, theta = theta)
  mu <- zsm(y=1:J, J = J, m = m, theta = theta)
  lpn <- log(sn) - log(sum(mu))
  if(log) return(lpn)
  else return(exp(lpn))
}

fun1 <- function(n, J, m, theta, x){
  gama <- m*(J-1)/(1-m)
  nu <- J+gama*(1-x)
  lambda <- gama*x
  ly <- exp(lchoose(J, n)+lgamma(n+lambda)-lgamma(lambda)+lgamma(nu-n)-lgamma(nu-J)+lgamma(lambda+nu-J)-lgamma(lambda+nu)+(theta-1)*log(1-x)-log(x))
  return(exp(ly))
}

zsm <- function(n, J, m, theta, x){
  f2 <- function(n, J, m, theta){
    integrate(f=fun1,lower=0,upper=1,n=n,J=J,m=m,theta=theta)$value
  }
  theta*sapply(n, f2, J, m, theta)
}

dzsm <- function(n, J, m, theta, log = FALSE){
  if (J <= 0) stop ("J must be great than zero")
  if (theta <= 0) stop ("theta must be great than zero")
  if (m > 1 | m < 0) stop ("m must be between zero and one")
  sn <- zsm(n, J = J, m = m, theta = theta)
  mu <- zsm(1:J, J = J, m = m, theta = theta)
  lpn <- log(sn)-log(sum(mu))
  if(log) return(lpn)
  else return(exp(lpn))
}
