dvolkov <- function(x, theta, m, J, log=FALSE, tol=1e-4, ...){
  X <- 1-m
  m.tilde <- J*m/(1-m)
  omega <- theta/m.tilde - log(1-X)
  f1 <- function(y,N){
    k <- log(theta)+(N*log(X)-lfactorial(N)) 
    f <- (lgamma(N+y)-lgamma(1+y))-omega*y 
    exp(k+f)
  }
  f2 <- function(ab){
   integrate(f1,0,Inf, N=ab, ...)$value
  }
  abunds <- 1
  vals <- c()
  vals[1] <- f2(abunds)
  abunds <- abunds+1
  vals[2] <- f2(abunds)
  while((sum(vals))/(sum(vals[-abunds]))>(1+tol)|abunds<max(x))
    {
      abunds <- abunds+1
      vals[abunds] <- f2(abunds)
      }
  ES <- theta*log((1-J*m)/(theta*(1-m))*log(m))
  if(ES>sum(vals)) Stot <- ES
  else Stot <- sum(vals)
  if(log)log(vals[x]/Stot)
  else vals[x]/Stot
}
