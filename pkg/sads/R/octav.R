octav <- function(x, oct, ...){
  if(is(x, "fitsad"))
    y <- x@data$x
  if(is(x,"fitrad"))
    y <- x@rad.tab$abund
  if(is(x,"numeric"))
    y <- x
  if(missing(oct)){
    oct <- 1:(ceiling(max(log2(y)))+1)
    if(any(y < 1)){
      octlower <- ceiling(min(log2((y)))+1):0
      oct <- c(octlower, oct)
    }
  }
  n <- 2^(oct-1)
  oc.class <- cut(y, breaks=c(0, n), labels=oct)
  res <- as.data.frame(table(oc.class))
  res$upper <- n
  names(res)[1] <- "octave"
  new("octav", res[,c(1,3,2)])
}
