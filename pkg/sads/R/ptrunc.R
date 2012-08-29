ptrunc <- function(f, q, trunc, ...){
  log.p <- FALSE
  lower.tail <- FALSE
  dots <- c(list(...))
  if("log.p"%in%names(dots)){
    log.p <- dots$log.p
    dots$log.p <- FALSE
  }
  if("log"%in%names(dots)){
    log.p <- dots$log
    dots$log <- FALSE
  }
  dots2 <- dots
  if("lower.tail"%in%names(dots2)){
    lower.tail <- dots2$lower.tail
    dots2$lower.tail <- FALSE
  }
  tt <- q
  pf <- get(paste("p", f, sep = ""), mode = "function")
  if(!missing(trunc)){
    aa <- rep(trunc, length(q))
    tt <- do.call(pf, c(list(q = apply(cbind(q, aa), 1, max)), dots))
    if(lower.tail){
      tt <- tt - do.call(pf, c(list(q = aa), dots))
    }
    tt <- tt/(1 - do.call(pf, c(list(q = aa), dots2)))
  }
  else{
    tt <- do.call(pf, c(list(q = q), dots))
  }
  if(log.p)return(log(tt)) else return(tt)
}
