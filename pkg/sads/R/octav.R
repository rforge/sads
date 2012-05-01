octav <- function(x, oct=1:(ceiling(max(log2((x))))+1), ...){
  n <- 2^(oct-1)
  oc.class <- cut(x,breaks=c(0,n), labels=oct)
  res <- as.data.frame(table(oc.class))
  res$upper <- n
  names(res)[1] <- "octave"
  new("octav",res[,c(1,3,2)])
}
