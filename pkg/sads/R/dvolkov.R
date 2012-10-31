dvolkov <- function(x, J, m, theta, log = FALSE){
  volkov <- function (J, m, theta){
    (gam <- m * (J - 1)/(1 - m))
    integrand <- function(y, n) {
      theta * exp(lgamma(J + 1) - lgamma(n + 1) - lgamma(J - 
        n + 1) + lgamma(gam) - lgamma(J + gam) + lgamma(n + 
        y) + lgamma(J - n + gam - y) - lgamma(1 + y) - lgamma(gam - 
        y) - y * theta/gam)
    }
    f <- function(n) {
      integrate(integrand, lower = 0, upper = gam , n = n)$value
    }
    out <- sapply(1:J, f)
    return(out)
  }
  Sj <- volkov(J = J, m = m, theta = theta)
  dsj <- Sj[x]/sum(Sj)
  if(log) return(log(dsj))
  else return(dsj)
}