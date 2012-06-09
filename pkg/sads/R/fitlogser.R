fitlogser <- function(x, size, rich, ...){
  if(!missing(x)){
    S <- length(x)
    N <- sum(x)
    if(!missing(size)|!missing(rich)){
      warning(paste("Model fitted with size = ",N," and rich = ",S," \n calculated from supplied abundances"))
    }
  }
  if(missing(x)&!missing(size)&!missing(rich)){
    S <- rich
    N <- size
  }
  if(missing(x)&missing(size)&missing(rich)){
    stop("Please provide size and species number or a vector of species abundances")
  }
  f1 <- function(a) {
        S + a * log((a/(a + N)))
    }
  sol <- uniroot(f1, interval = c(1/N, N))
  alfa=sol$root
  X <- N/(N+alfa)
  if(!missing(x)){
    LL <- function(alpha)-sum(dls(x,N,alpha,log=T),...)
    result <- mle2(LL,start=list(alpha=alfa),data=list(x=x))
    new("fitsad",result, sad="ls")
  }
  else new("fitsad",coef=c(alpha=alfa),fullcoef=c(alpha=alfa),sad="ls")
}
