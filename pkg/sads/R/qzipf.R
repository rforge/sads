qzipf<-function(p, N, s){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- pzipf(U2, N, s)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(pzipf(min(a1, a2), N, s) < U1 & U1 <= pzipf(max(a1, a2), N, s)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=N)) #problemas no max
    if(U1 <= pzipf(0, N, s)){
      d[i] <- 0
    } else if(U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}