fitpoilog <- function(x, trunc=0,...){
  pl.par <- poilogMLE(x,startVals=c(mu=mean(log(x))+log(0.5),sig=sd(log(x))))$par
  LL <- function(mu,sig) -sum(dpoilog(x,mu,sig,trunc=trunc,log=TRUE))
  new("fitsad",mle2(LL,start=as.list(pl.par),data=list(x=x),...),sad="poilog")
}
