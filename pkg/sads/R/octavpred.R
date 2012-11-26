octavpred <- function(object, sad, rad, coef, trunc, oct, S, N, ...){
  dots <- list(...)
  if(missing(trunc)) trunc <- NA
  if(missing(object)){
    if((missing(sad)&&missing(rad))||missing(coef)||missing(oct)||missing(S)||missing(N))
      stop("please provide 'object' or input  'sad' or 'rad' plus 'coef', 'oct', 'S' , and 'N' ")
  }
  else{
    if(class(object)=="fitsad"||class(object)=="fitrad"){
      coef <- as.list(bbmle::coef(object))
      trunc <- object@trunc
      if(class(object)=="fitsad"){
        sad <- object@sad
        x <- object@data$x
      }
      if(class(object)=="fitrad"){
        rad <- object@rad
        x <- object@rad.tab$abund
      }
    }
    else if(class(object)=="numeric"){
      if((missing(sad)&&missing(rad))||missing(coef))
        stop("please input 'sad' or 'rad' and 'coef' ")
      else x <- object
    }
    S <- length(x)
    N <- sum(x)
    if(missing(oct)){
      oct <- 1:(ceiling(max(log2(x)))+1)
      if(any(x < 1)){
        octlower <- ceiling(min(log2((x)))+1):0
        oct <- c(octlower, oct)
      }
    }
  }
  n <- 2^(oct-1)
  if(!missing(sad)){
    if(!is.na(trunc)){
      if(sad == "ls")
        Y <- do.call(ptrunc, c(list(f=sad, q = n, coef=c(list(N = N),coef),trunc = trunc),dots))
      else if(sad == "mzsm"||sad=="volkov")
        Y <- do.call(ptrunc, c(list(f=sad, q = n, coef = c(list(J = N),coef), trunc = trunc),dots))
      else
        Y <- do.call(ptrunc, c(list(sad, q = n, coef = coef, trunc = trunc), dots))
    }
    else{
      psad <- get(paste("p",sad,sep=""),mode="function")
      if(sad == "ls")
        Y <- do.call(psad, c(list(q = n, N = N),coef,dots))
      else if(sad == "mzsm" || sad == "volkov")
        Y <- do.call(psad, c(list(q = n, J = N),coef,dots))
      else
        Y <- do.call(psad, c(list(q = n),coef,dots))
    }
    Y <- c(Y[1], diff(Y))*S
  }
  ## for rads:
  else if(!missing(rad)){
    if(!is.na(trunc)){
      ab <- do.call(dtrunc, c(list(f=rad, q = 1:S, coef=coef,trunc = trunc),dots))*N
    }
    else{
      drad <- get(paste("d",rad,sep=""),mode="function")
      ab <- do.call(drad, c(list(x=1:S),coef,dots))*N
    }
    Y <- as.numeric(table(cut(ab,breaks=c(2^(min(oct)-2),n))))
  }  
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}

