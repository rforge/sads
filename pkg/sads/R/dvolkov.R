dvolkov <- function(x, theta, m, J, log=FALSE){
  vals <- volkov(J,c(theta,m))
  Stot <- sum(vals)
  if(log)log(vals[x]/Stot)
  else vals[x]/Stot
}
