radpred <- function(object, sad, rad, coef, trunc, distr, S, N, ...){
  dots <- list(...)
  if(missing(trunc)) trunc <- NA
  if(missing(object)){
    if((missing(sad)&&missing(rad))||missing(coef)||missing(distr)||missing(S)||missing(N))
      stop("please provide 'object' or input  'sad' or 'rad' plus 'coef', 'distr', 'S' , and 'N' ")
  }
  else{
    if(class(object)=="fitsad"||class(x=object)=="fitrad"){
      coef <- as.list(bbmle::coef(object))
      trunc <- object@trunc
      distr <- object@distr
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
  }
  ## For rad models
  if(!missing(rad)){
    y <- 1:S
    if(!is.na(trunc)){
      ab <- do.call(dtrunc, c(list(rad, x = y, coef = coef, trunc = trunc), dots))*N
    }
    else{
      drad <- get(paste("d", rad, sep=""),  mode = "function")
      ab <- do.call(drad, c(list(x = y), coef, dots))*N
    }
  }
### For sad models    
  else if(!missing(sad)){
    if (distr == "D"){
      y <- 1:N
      if(!is.na(trunc)){
        if(sad=="ls") 
          X <- do.call(ptrunc, c(list(sad, q = y, coef = c(list(N = N),coef),
                                      lower.tail=F, trunc = trunc)))
        else if(sad=="mzsm"||sad=="volkov") 
          X <- do.call(ptrunc, list(sad, q = y, coef = c(list(J = N), coef),
                                    lower.tail=F, trunc = trunc))
        else
          X <- do.call(ptrunc, list(sad, q = y, coef = coef, lower.tail=F, trunc = trunc))
      }
      else {
        psad <- get(paste("p", sad, sep=""), mode = "function")
        if(sad=="ls")
          X <- do.call(psad, c(list(q = y, lower.tail = F, N = N),coef))
        else if(sad=="mzsm"||sad=="volkov")
          X <- do.call(psad, c(list(q = y, lower.tail = F, J = N), coef))
        else
          X <- do.call(psad, c(list(q = y, lower.tail = F), coef))
      }
      f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
      ab <- f1(ppoints(S))
    }
    else if(distr == "C"){
      Y <- ppoints(S)
      if(!is.na(trunc)){
        ab <- do.call(qtrunc, list(sad, p = Y, coef = coef, lower.tail=F, trunc = trunc))
      }
      else{
        qsad <- get(paste("q", sad, sep=""), mode = "function")
        ab <- do.call(qsad, c(list(p = Y, lower.tail = F), coef))
      }
    }
    else
      stop("unsupported distribution")
  }
  new("rad", data.frame(rank=1:S, abund=ab))
}
