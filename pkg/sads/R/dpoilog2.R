dpoilog2 <- function(x, mu, sig, trunc=0, log = FALSE){
  if(!is.null(trunc)){
    if(any(x<=trunc)) stop("at least one x larger than truncation value")
    iconst <- 1-sum(poilog::dpoilog(0:trunc, mu=mu, sig=sig))
  }
  else iconst <- 1
  y <- (poilog::dpoilog(x, mu, sig))/iconst
  if(log) return(log(y))
  else return(y)
}