qqsad <- function(object, ...){
  dots <- list(...)
  rank <- sort(object@data$x)
  S <- length(object@data$x)
  p <- ppoints(S)
  if(!is.na(object@trunc)){
    q <- do.call(qtrunc, c(list(object@sad, p = p, coef = as.list(object@coef), trunc = object@trunc), dots))
  }else{
    qsad <- get(paste("q", object@sad, sep=""), mode = "function")
    q <- do.call(qsad, c(list(p = p), as.list(object@coef), dots))
  }
  plot(q, rank, main = "Q-Q plot", xlab="Theoretical Quantile", ylab="Sample Quantiles")
  abline(0, 1, col = "red", lty = 2)
}