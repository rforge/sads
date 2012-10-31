system("R CMD SHLIB zsm.c")
dyn.load("zsm.so")


# ZSM ---------------------------------------------------------------------
zsm <- function(n, J, m, theta, precisao = 0.01){
  lengn <- length(n)
  res <- rep(0, lengn)
  fun <- .C("sn", as.integer(n), as.integer(lengn), as.integer(J), as.double(m),
            as.double(theta), as.double(precisao), result = res)
  fun[["result"]]
}

integr2 <- function(n, J, m, theta, precisao = 0.01){
  res <- 0
  fun <- .C("Intgrl1", as.integer(n), as.integer(J), as.double(m),
            as.double(theta), as.double(precisao), result = res)
  fun[["result"]]
}

zsm2 <- function(n, J, m, theta, ...){
  f1 <- function(n){
    integr2(n, J, m, theta,...)
  }
  theta*sapply(n, f1)
}

system.time((teste <- alonso10(1:150, J=150, m=0.05, theta=2)))
sum(teste)

system.time((teste2 <- zsm2(1:150, J=150, m=0.05, theta=2)))
sum(teste2)

system.time((teste3 <- zsm(1:150, 150, 0.05, 2)))
sum(teste3)

plot(teste, teste2)
abline(0,1)
plot(1:150, teste)

alonso14 <- function(n, J, m, theta, log = FALSE,...){
  all.values <- alonso10(1:J, J, m, theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}

dzsm <- function(n, J, m, theta, log = FALSE,...){
  all.values <- zsm2(1:J, J, m, theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}

system.time((teste3 <- alonso14(1:150, 150, 0.05, 2)))
sum(teste3)
system.time((teste4 <- dzsm(1:150, 150, 0.05, 2)))
sum(teste4)

plot(teste3, teste4)
abline(0,1)
plot(1:150, teste3)

(alonso14(n=1:5, J=15000, m=0.77, theta=41))
(dzsm(n=1:5, J=15000, m=0.77, theta=41))


# MZSM --------------------------------------------------------------------
mzsm <- function(n, J, theta){
  lengn <- length(n)
  res <- rep(0, lengn)
  fun <- .C("msn", as.integer(n), as.integer(lengn), as.integer(J), as.double(theta), result = res)
  fun[["result"]]
}

system.time((teste <- mzsm(1:150, 150, 0.9)))
sum(teste)
head(teste)

####theta < 1 e J > 100
system.time((teste2 <- alonso10(1:150, 150, 0.999, 0.9)))
sum(teste2)
tail(teste) 

plot(teste, teste2)
abline(0,1)
plot(1:150, teste)

dmzsm <- function(n, J, theta, log = FALSE,...){
  all.values <- mzsm(1:J, J, theta)
  all.values[all.values == Inf] <- NaN
  lprobs <- log(all.values[n])-log(sum(all.values, na.rm=TRUE))
  if(log) lprobs
  else exp(lprobs)
}

system.time((teste3 <- alonso14(1:150, 150, 0.99, 0.9)))
sum(teste3)
system.time((teste4 <- dmzsm(1:150, 150, 0.9)))
sum(teste4, na.rm=T)

plot(teste3, teste4)
abline(0,1)
plot(1:150, teste4, type ="l")
lines(1:150, teste3, col = 2)