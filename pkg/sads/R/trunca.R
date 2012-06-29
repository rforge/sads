trunca <- function(f, x, trunc,...){
  dots <- list(...)
  p <- get(paste("p", deparse(substitute(f)), sep=""), mode = "function")
  d <- get(paste("d", deparse(substitute(f)), sep=""), mode = "function")
  tt <- rep(0, length(x))
  if (!missing(trunc)){
    tt[x>= trunc] <- do.call(d, c(list(x=x[x>=trunc]),dots))/(1-do.call(p,c(list(q=trunc),dots)))    
  }
  else{
        tt <- do.call(d, c(list(x=x),dots))        
  }
  return(tt)
}
