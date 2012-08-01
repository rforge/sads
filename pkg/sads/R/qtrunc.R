qtrunc <- function(f, p, trunc, ...){
  dots <- c(list(...))
  tt <- p
  pf <- get(paste("p", f, sep = ""), mode = "function")
  qf <- get(paste("q", f, sep = ""), mode = "function")
  if(!missing(trunc)){
    aa <- do.call(pf, c(list(q = trunc), dots)) + p*(1 - do.call(pf, c(list(q = trunc), dots)))
    tt <- do.call(qf, c(list(p = aa), dots))
  } else{
    tt <- do.call(qf, c(list(p = p), dots))
  }
  return(tt)
}