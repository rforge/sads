pred.logser=function(x,alpha,size,rich){
  if(is.null(alpha)&is.null(size)&is.null(rich)){
    stop("Please provide at least two of these: alpha, size, rich")
  }
  if(missing(alpha)){
    alpha <- as.numeric(coef(fitlogser(size=size,rich=rich)))
  }
  if(missing(size)){
    size <- alpha*exp(rich/alpha) - alpha
  }
  X <- size/(alpha + size)
  alpha*(X^x)/x
}
