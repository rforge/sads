octav <- function(x, oct, ...){
  if(is(x, "fitsad")){
    if(missing(oct)){
      oct <- 1:(ceiling(max(log2(x@data$x)))+1)
      if(any(x@data$x < 1)){
        octlower <- ceiling(min(log2((x$data$x)))+1):0
        oct <- c(octlower, oct)
      }
    }
    n <- 2^(oct-1)
    oc.class <- cut(x@data$x, breaks=c(0, n), labels=oct)
    res <- as.data.frame(table(oc.class))
    res$upper <- n
    names(res)[1] <- "octave"
  }else{
    if(missing(oct)){
      oct <- 1:(ceiling(max(log2(x)))+1)
      if(any(x < 1)){
        octlower <- ceiling(min(log2((x)))+1):0
        oct <- c(octlower, oct)
      }
    }
    n <- 2^(oct-1)
    oc.class <- cut(x, breaks=c(0, n), labels=oct)
    res <- as.data.frame(table(oc.class))
    res$upper <- n
    names(res)[1] <- "octave"
  }
  new("octav", res[,c(1,3,2)])
}
