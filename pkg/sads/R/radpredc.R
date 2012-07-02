radpredc <- function(x, sad, coef,...){
  if(missing(coef))dots <- list(...)
  else dots <- c(coef, list(...)
  qsad <- get(paste("q", deparse(substitute(sad)), sep=""), mode = "function")
  S <- length(x)
  Y <- ppoints(S)
  ab <- do.call(qsad,c(list(p=Y,lower.tail=F),dots))
  new("rad",data.frame(rank=1:S, abund=ab))
}
