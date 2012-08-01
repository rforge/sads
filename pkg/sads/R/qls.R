qls<-function(p, N, alpha){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- pls(U2, N, alpha)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(pls(min(a1, a2), N, alpha) < U1 & U1 <= pls(max(a1, a2), N, alpha)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=N))
    if(U1 <= pls(1, N, alpha)){
      d[i] <- 1
    } else if (U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}