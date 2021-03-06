dsad <- function(y,a,sad,log=FALSE,upper=0.9999,...){
  qcom <- paste("q",deparse(substitute(sad)),sep="")
  dcom <- paste("d",deparse(substitute(sad)),sep="")
  dots <- c(as.name("n"),list(...))
  uplim <- do.call(qcom,c(upper,list(...)))
  f1 <- function(z){
    f2 <- function(n){
      t1 <- do.call(dcom,dots)
      t2 <- dpois(z,a*n)
      ifelse(t1==0|t2==0,0,t1*t2*1e12)
    }
  integrate(f2,0,uplim,rel.tol=sqrt(.Machine$double.eps),subdivisions=500)$value
  }
  res <- sapply(y,f1)/(1e12*upper)
  if(log)log(res) else res
}

