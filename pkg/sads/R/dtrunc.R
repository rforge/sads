dtrunc <- function(f, x, trunc, log = FALSE, ...){
  dots <- c(list(...))
  pf <- get(paste("p", f, sep=""), mode = "function")
  df <- get(paste("d", f, sep=""), mode = "function")
  tt <- rep(0, length(x))
  if (!missing(trunc)){
    tt[x > trunc] <- do.call(df, c(list(x = x[x>trunc]), dots))/(1 - do.call(pf, c(list(q = trunc), dots)))    
  } else{
    tt <- do.call(df, c(list(x = x), dots))
  }
  if (log) tt <- log(tt)
  return(tt)
}