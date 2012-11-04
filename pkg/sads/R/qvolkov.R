qvolkov<-function(p, theta, m, J, lower.tail = TRUE, log.p = FALSE, cor.lin=FALSE,...){
  if (length(theta) > 1 | length(m) > 1 |length(J) > 1) stop("vectorization of theta, m and J is not implemented")
  if (!all(is.finite(c(J, theta, m)))) stop("all parameters should be finite")
  if (J <= 0)  stop("J must be larger than zero")
  if (theta <= 0) stop("theta must be larger than zero")
  if(m < 0 | m > 1) stop("m must be between zero and one")
  if (log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  Y <- 1:J
  X <- pvolkov(Y, theta=theta, m=m, J=J,...)
  if(cor.lin){
    k <- (1-Y[J])/sum(Y==Y[J])
    Y[Y==Y[J]] <- cumsum(rep(k,sum(Y==Y[J])))+Y[J]
  }
  f1 <- approxfun(X,Y,method="constant",f=0,yleft=1, yright=J)
  f1(p)
}
