ptrunc <- function(f, q, trunc, ...){
  dots <- c(list(...))
  tt <- q
  pf <- get(paste("p", f, sep = ""), mode = "function")
  if(!missing(trunc)){
    aa <- rep(trunc, length(q))
    tt <- do.call(pf, c(list(q = apply(cbind(q, aa), 1, max)), dots))
    tt <- tt - do.call(pf, c(list(q = aa), dots))
    tt <- tt/(1 - do.call(pf, c(list(q = aa), dots)))
  } else{
    tt <- do.call(pf, c(list(q = q), dots))
  }
  return(tt)
}