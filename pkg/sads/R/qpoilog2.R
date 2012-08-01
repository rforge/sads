qpoilog2<-function(p, mu, sig, S=30){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- ppoilog2(U2, mu, sig)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(ppoilog2(min(a1, a2), mu, sig) < U1 & U1 <= ppoilog2(max(a1, a2), mu, sig)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=S))
    if(U1 <= ppoilog2(0, mu, sig)){
      d[i] <- 0
    } else if(U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}