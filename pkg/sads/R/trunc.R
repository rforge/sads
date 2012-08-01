trunc<-function(f, x, trunc, log = FALSE, ...){
  #dots <-c(list(...))
  first<-strsplit(f, NULL)[[1]][1]
  fun<- paste(strsplit(f, NULL)[[1]][-1], collapse="")
  switch(first, 
         d = dtrunc(fun, x, trunc, log, ...), 
         p = ptrunc(fun, x, trunc, ...),
         q = qtrunc(fun, x, trunc, ...))
}

#teste d
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

#teste p
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

#teste q
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