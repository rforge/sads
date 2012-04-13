pred.logser=function(x,alpha=NULL,size=NULL,rich=NULL){
  if(is.null(alpha)&is.null(size)&is.null(rich)){
    stop("Please provide at least two of these: alpha, size, rich")
  }
  if(is.null(alpha)){
    alpha <- fit.logser(size=size,rich=rich)$summary["alpha"]
  }
  if(is.null(size)){
    size <- alpha*exp(rich/alpha) - alpha
  }
  X <- size/(alpha + size)
  alpha*(X^x)/x
}
