octavpred <- function(object, sad, coef, trunc, oct, S, N, ...){
  dots <- list(...)
  if(missing(trunc)) trunc <- NA
  if(missing(object)){
    if(missing(S)||missing(N)||missing(sad)||missing(coef)||missing(oct))
      stop("please provide 'object' or input 'sad', 'coef', 'oct', 'S' , and 'N' ")
  }
  else{
    if(class(object)=="fitsad"||class(x)=="fitrad"){
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
    else if(class(object)=="numeric")
      x <- object
    S <- length(x)
    N <- sum(x)
    if(missing(oct)){
      oct <- 1:(ceiling(max(log2(x)))+1)
      if(any(x < 1)){
        octlower <- ceiling(min(log2((x)))+1):0
        oct <- c(octlower, oct)
      }
    }
    else if(!missing(oct)){
      oct <- 1:oct
    }
    else stop("input a fitsad, fitsad or abundance vector or input oct")
    n <- 2^(oct-1)
    if(!missing(sad)){
      if(!is.na(trunc)){
        if(sad == "ls")
          Y <- do.call(ptrunc, c(list(f=sad, q = n, coef=c(list(N = N),coef),trunc = trunc),dots))
        else if(sad == "mzsm"||"volkov")
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
    ##Incluir os comandos para rads:
    ##else if(!missing(rad)){
    ##}
  }  
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}

