dpoilog2 <- function(x, mu, sig, log = FALSE){
  y <- poilog::dpoilog(x, mu, sig)
  if(log) return(log(y))
  else return(y)
}