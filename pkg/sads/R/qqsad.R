qqsad <- function(object, ...){
  rank <- sort(object@data$x)
  S <- length(object@data$x)
  p <- ppoints(S)
  qsad <- get(paste("q", object@sad, sep=""), mode = "function")
  q <- qsad(p)
  plot(rank, q, col = "red")
  abline(0, 1)
}