library(sads); library(VGAM) #distrb. zipf (dzipf, pzipf)
source("trunc.R")
source("qls.R")
source("qpoilog.R")
source("qzipf.R")

source("dpoilog2.R")
source("ppoilog2.R")
source("qpoilog2.R")
#a diferenca das funcoes poilog2 para as funcoes poilog e a opcao de truncagem

# teste trunc, dtrunc, ptrunc e qtrunc
###############################
# Distrb continua - Normal    #
###############################
x<- -4:4
dnorm(x)
dtrunc("norm", x)
trunc("dnorm", x)

x<- -2:8
dnorm(x, 4, 2)
dtrunc("norm", x, mean =4, sd = 2)
trunc("dnorm", x, mean= 4, sd = 2)

dtrunc("norm", x, mean =4, sd = 2, trunc=0)
trunc("dnorm", x, mean= 4, sd = 2, trunc=0)

q<--4:4
pnorm(q)
ptrunc("norm", q)
trunc("pnorm", q)

q<- -2:8
pnorm(q, 4, 2)
ptrunc("norm", q, mean =4, sd = 2)
trunc("pnorm", q, mean= 4, sd = 2)

ptrunc("norm", q, mean =4, sd = 2, trunc=0)
trunc("pnorm", q, mean= 4, sd = 2, trunc=0)
#undebug(ptrunc)

p<-seq(0, 1, 0.01)
qnorm(p)
qtrunc("norm", p)
trunc("qnorm", p)

qnorm(p, 4, 2)
qtrunc("norm", p, mean =4, sd = 2)
trunc("qnorm", p, mean= 4, sd = 2)

qtrunc("norm", p, mean =4, sd = 2, trunc=0)
trunc("qnorm", p, mean= 4, sd = 2, trunc=0)

#################################
# Distr. discreta - Poisson    ##
#################################
x<- 0:10
dpois(x, 2)
dtrunc("pois", x, lambda = 2)
trunc("dpois", x, lambda = 2)

dtrunc("pois", x, lambda = 2, trunc = 0)
trunc("dpois", x, lambda = 2, trunc = 0)

q<-0:10
ppois(q, lambda = 2)
p1<-ptrunc("pois", q, lambda = 2)
trunc("ppois", q, lambda = 2)

ptrunc("pois", q, lambda = 2, trunc = 2)
p1<-trunc("ppois", q, lambda = 2, trunc = 2)

p<-seq(0, 1, 0.01)
q<- qpois(p, lambda = 2);q
q<- qtrunc("pois", p, lambda = 2);q
q<- trunc("qpois", p, lambda = 2);q

q<- qtrunc("pois", p, lambda = 2, trunc=2)
q<- trunc("qpois", p, lambda = 2, trunc=2)

###################################
# Distr. discreta - Logserie      #
###################################
x<- 1:20
dls(x, N = 10, alpha = 2)
dtrunc("ls", x, N = 10, alpha = 2)
trunc("dls", x, N = 10, alpha = 2)

dtrunc("ls", x, N = 10, alpha = 2, trunc = 2)
trunc("dls", x, N = 10, alpha = 2, trunc = 2)
#debug(dtrunc)

q<-1:20
pls(q, N = 10, alpha = 2)
ptrunc("ls", q, N = 10, alpha = 2)
trunc("pls", q, N = 10, alpha = 2)

ptrunc("ls", q, N = 10, alpha = 2, trunc = 2)
trunc("pls", q, N = 10, alpha = 2, trunc = 2)
#undebug(ptrunc)

p<- seq(0, 0.99, 0.01)
qls(p, N = 10, alpha= 2)
qtrunc("ls", p, N = 10, alpha = 2)
trunc("qls", p, N = 10, alpha = 2)

qtrunc("ls", p, N = 10, alpha = 2, trunc = 2)
trunc("qls", p, N = 10, alpha = 2, trunc = 2)

p<- seq(0, 0.999, 0.001)
q<- qls(p, N = 10, alpha= 2)
p1 <- pls(q, N = 10, alpha = 2)
plot(q, p1)
q<- qtrunc("ls", p, N = 10, alpha = 2)
p1<- ptrunc("ls", q, N = 10, alpha = 2)
points(q, p1, col=2, pch="#")
q<- trunc("qls", p, N = 10, alpha = 2)
p1<- trunc("pls", q, N = 10, alpha = 2)
points(q, p1, col=3, pch="*")

q<- qtrunc("ls", p, N = 10, alpha = 2, trunc = 2)
p1<- ptrunc("ls", q, N = 10, alpha = 2, trunc = 2)
plot(q, p1)
q<- trunc("qls", p, N = 10, alpha = 2, trunc = 2)
p1<- trunc("pls", q, N = 10, alpha = 2, trunc = 2)
points(q, p1, col=2, pch="*")

#########################################
#Distr. discreta - Poisson lognormal   ##
#########################################
#Obs.: A funcao dpoilog tem que ser truncada em zero? (Funcao dpoilog no sads e truncada em zero)
x<- 1:20
dpoilog(x, mu = 1, sig = 1)
dtrunc("poilog", x, mu = 1, sig = 1)
trunc("dpoilog", x, mu = 1, sig = 1)

dtrunc("poilog", x, mu = 1, sig = 1, trunc = 2)
trunc("dpoilog", x, mu = 1, sig = 1, trunc = 2)

q<-1:20
ppoilog(q, mu = 1, sig = 1)
ptrunc("poilog", q, mu = 1, sig = 1)
trunc("ppoilog", q, mu = 1, sig = 1)

ptrunc("poilog", q, mu = 1, sig = 1, trunc = 2)
trunc("ppoilog", q, mu = 1, sig = 1, trunc = 2)
#undebug(ptrunc)

p<- seq(0, 0.99, 0.01)
qpoilog(p, mu = 1, sig= 1)
qtrunc("poilog", p, mu = 1, sig = 1)
trunc("qpoilog", p, mu = 1, sig = 1)

qtrunc("poilog", p, mu = 1, sig = 1, trunc = 2)
trunc("qpoilog", p, mu = 1, sig = 1, trunc = 2)

p<- seq(0, 0.99, 0.01)
q<-qpoilog(p, mu = 1, sig= 1)
p1<-ppoilog(q, mu = 1, sig = 1)
plot(q, p1)
q<- qtrunc("poilog", p, mu = 1, sig = 1)
p1<- ptrunc("poilog", q, mu = 1, sig = 1)
points(q, p1, col=2, pch="#")
q<- trunc("qpoilog", p, mu = 1, sig = 1)
p1<- trunc("ppoilog", q, mu = 1, sig = 1)
points(q, p1, col=3, pch="*")

q<-qtrunc("poilog", p, mu = 1, sig = 1, trunc = 2)
p1<-ptrunc("poilog", q, mu = 1, sig = 1, trunc = 2)
plot(q, p1)
q<-trunc("qpoilog", p, mu = 1, sig = 1, trunc = 2)
p1<-trunc("ppoilog", q, mu = 1, sig = 1, trunc = 2)
points(q, p1, col=2, pch="*")

### ppoilog e dpoilog sem opcao de truncagem
x<- 0:20
dpoilog2(x, mu = 1, sig = 1)
dtrunc("poilog2", x, mu = 1, sig = 1)
trunc("dpoilog2", x, mu = 1, sig = 1)

dtrunc("poilog2", x, mu = 1, sig = 1, trunc = 0)
trunc("dpoilog2", x, mu = 1, sig = 1, trunc = 0)
#undebug(dtrunc)

q<-0:20
ppoilog2(q, mu = 1, sig = 1)
ptrunc("poilog2", q, mu = 1, sig = 1)
trunc("ppoilog2", q, mu = 1, sig = 1)

ptrunc("poilog2", q, mu = 1, sig = 1, trunc = 0)
trunc("ppoilog2", q, mu = 1, sig = 1, trunc = 0)
#undebug(ptrunc)

p<- seq(0, 0.99, 0.01)
qpoilog2(p, mu = 1, sig= 1)
qtrunc("poilog2", p, mu = 1, sig = 1)
trunc("qpoilog2", p, mu = 1, sig = 1)

qtrunc("poilog2", p, mu = 1, sig = 1, trunc = 2)
trunc("qpoilog2", p, mu = 1, sig = 1, trunc = 2)

p<- seq(0, 0.99, 0.01)
q<-qpoilog2(p, mu = 1, sig= 1)
p1<-ppoilog2(q, mu = 1, sig = 1)
plot(q, p1)
q<- qtrunc("poilog2", p, mu = 1, sig = 1)
p1<- ptrunc("poilog2", q, mu = 1, sig = 1)
points(q, p1, col=2, pch="#")
q<- trunc("qpoilog2", p, mu = 1, sig = 1)
p1<- trunc("ppoilog2", q, mu = 1, sig = 1)
points(q, p1, col=3, pch="*")

q<-qtrunc("poilog2", p, mu = 1, sig = 1, trunc = 2)
p1<-ptrunc("poilog2", q, mu = 1, sig = 1, trunc = 2)
plot(q, p1)
q<-trunc("qpoilog2", p, mu = 1, sig = 1, trunc = 2)
p1<-trunc("ppoilog2", q, mu = 1, sig = 1, trunc = 2)
points(q, p1, col=2, pch="*")

########################
# Distribuicao ZIPF   ##
########################
x<- 0:15
dzipf(x, N = 30, s = 2)
dtrunc("zipf", x, N = 30, s = 2)
trunc("dzipf", x, N = 30, s = 2)

dtrunc("zipf", x, N = 30, s = 2, trunc = 1)
trunc("dzipf", x, N = 30, s = 2, trunc = 1)

q<-0:15
pzipf(q, N = 30, s = 2)
ptrunc("zipf", q, N = 30, s = 2)
trunc("pzipf", q, N = 30, s = 2)

ptrunc("zipf", q, N = 30, s = 2, trunc = 1)
trunc("pzipf", q, N = 30, s = 2, trunc = 1)

p<-seq(0, 1, 0.01)
q<- qzipf(p, N = 30, s = 2); q
p1<-pzipf(q[-101], N = 30, s = 2); p1
plot(q[-101], p1)
q<- qtrunc("zipf", p, N = 30, s = 2);q
p1<-ptrunc("zipf", q[-101], N = 30, s = 2)
points(q[-101], p1, col=2, pch="#")
q<- trunc("qzipf", p, N = 30, s = 2);q
p1<- trunc("pzipf", q[-101], N = 30, s = 2)
points(q[-101], p1, col=3, pch="*")

q<- qtrunc("zipf", p, N = 30, s = 2, trunc=1)
p1<-ptrunc("zipf", q[-101], N = 30, s = 2, trunc = 1)
plot(q[-101], p1)
q<- trunc("qzipf", p, N = 30, s = 2, trunc=1)
p1<-trunc("pzipf", q[-101], N = 30, s = 2, trunc = 1)
points(q[-101], p1, col=2, pch="*")

########################################################################
### Verificacao das funcoes ppoilog e pls                          #####
########################################################################
ppoilog3 <- function(q, mu, sig, lower.tail=TRUE, log.p=FALSE){
  y <- NULL
  y[1] <- dpoilog(q[1], mu, sig)
  if(length(q)>1){
    for(i in 2:length(q)){
      y[i] <- y[i-1] + dpoilog(q[i], mu, sig)
    }
  }
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}

q<-1:1000
system.time(d1<-ppoilog(q, 10, 13))  # media de 5.464
system.time(d2<-ppoilog2(q, 10, 13)) # media de 5.397
system.time(d3<-ppoilog3(q, 10, 13)) # media de 0.188

q<-1:10
trunc("ppoilog", q, mu=10, sig=13, trunc=2)
trunc("ppoilog2", q, mu=10, sig=13, trunc=2)
trunc("ppoilog3", q, mu=10, sig=13, trunc=2)

#A funcao ppoilog3 gasta menos tempo entretanto nao serve para a funcao trunc


pls2 <- function(q, N, alpha, lower.tail=TRUE, log.p=FALSE){
  y <- NULL
  y[1] <- dls(q[1], N, alpha)
  for(i in 2:length(q)){
    y[i] <- y[i-1] + dls(q[i], N, alpha)
  }
  if(!lower.tail) y <- 1-y
  if(log.p) y <- log(y)
  return(y)
}

q<-1:1000
system.time(d1<-pls(q, 1000, 1.3))  # Media de 0.168
system.time(d2<-pls2(q, 1000, 1.3)) # Media de 0.032

# mesmo problema da ppoilog3 nao funciona da trunc


########################################################################
#### Teste para a funcao de busca quantilica                   #########
########################################################################
### Poisson
qpois2<-function(p, lambda){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- ppois(U2, lambda)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(ppois(min(a1, a2), lambda) < U1 & U1 <= ppois(max(a1, a2), lambda)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=3*lambda)) #problemas no max
    if(U1 <= ppois(0, lambda)){
      d[i] <- 0
    } else if (U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}

p<-seq(0, 1, 0.01)
plot(p, qpois(p, lambda=2))
points(p, qpois2(p, lambda=2), col=2, pch="*")

q<-0:10
pv<-ppois(q, 4)
qpois(pv, 4)
qpois2(pv, 4)
#undebug(qpois2)

p<-seq(0, 1, 0.001)
q1<-qpois(p, 4)
q2<-qpois2(p, 4)
plot(p, q1)
points(p, q2, col=2)

q<-0:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)

q<-3:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)

q<-0:30
p<-ppois(q, 4)
plot(q, qpois2(p, 4), col=2)
points(q, qpois(p, 4))

### Geometrica
qgeom2<-function(p, prob){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- pgeom(U2, prob)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(pgeom(min(a1, a2), prob) < U1 & U1 <= pgeom(max(a1, a2), prob)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=20)) 
    if(U1 <= pgeom(0, prob)){
      d[i] <- 0
    } else if (U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}

p<-seq(0, 1, 0.01)
plot(p, qgeom(p, prob=1/2))
points(p, qgeom2(p, prob=1/2), col=2, pch = "*")

q<-0:10
pv<-pgeom(q, 1/4)
qgeom(pv, 1/4)
qgeom2(pv, 1/4)

p<-seq(0, 1, 0.001)
q1<-qgeom(p, 1/4)
q2<-qgeom2(p, 1/4)
plot(p, q1)
points(p, q2, col=2)

q<-0:15
p<-pgeom(q, 1/4)
qgeom(p, 1/4)
qgeom2(p, 1/4)

q<-3:15
p<-pgeom(q, 1/4)
qgeom(p, 1/4)
qgeom2(p, 1/4)

q<-0:130
p<-pgeom(q, 1/4)
plot(q, qgeom(p, 1/4), type="l")
lines(q, qgeom2(p, 1/4), col=2)

###############################################################################
#approxfun para funcao quantilica                                           ###
###############################################################################
#Poisson
qpois2<-function(p, lambda, mx=100){
  q  <- 0:mx
  tt <- p
  pv <- ppois(q, lambda)
  v <- approxfun(pv, q, method="linear") #linear e default
  for(i in 1:length(tt)){
    if(tt[i] <= ppois(0, lambda)){
      tt[i] <- 0
    } else if(tt[i] >= 0.999999999999999){
      tt[i] <- Inf
    }else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.001)
q1<-qpois(p, 4)
q2<-qpois2(p, 4)
plot(p, q1)
points(p, q2, col=2)

plot(q, pv, main = "approx(.) and approxfun(.)")
points(approx(q, pv), col = 2, pch = "*")
points(approx(q, pv, method = "constant"), col = 4, pch = "*")
points(q, pv, col = 3)
points(approx(q, pv, xout = 3, method = "constant"), col = 2, pch = "#")

q<-0:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)

q<-3:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)

q<-0:150
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4, mx=150)

#Geometrica
qgeom2<-function(p, prob, mx=100){
  q  <- 0:mx
  tt <- p
  pv <- pgeom(q, prob)
  v <- approxfun(pv, q, method="linear")
  for(i in 1:length(tt)){
    if(tt[i] <= pgeom(0, prob)){
      tt[i] <- 0
    } else if(tt[i] >= 0.999999999999999){
      tt[i] <- Inf
    }else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.001)
q1<-qgeom(p, 1/4)
q2<-qgeom2(p, 1/4)
plot(p, q1)
points(p, q2, col=2)

q<-0:15
p<-pgeom(q, 1/4)
qgeom(p, 1/4)
qgeom2(p, 1/4)

q<-3:15
p<-pgeom(q, 1/4)
qgeom(p, 1/4)
qgeom2(p, 1/4)

q<-0:150
p<-pgeom(q, 1/4)
qgeom(p, 1/4)
qgeom2(p, 1/4, mx=150)


###########################################################################
### Teste das funcoes busca e de interpolacao                         #####
###########################################################################
### Poisson
# Busca
qpois2<-function(p, lambda){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- ppois(U2, lambda)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(ppois(min(a1, a2), lambda) < U1 & U1 <= ppois(max(a1, a2), lambda)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=lambda/2, max=lambda))
    if(U1 <= ppois(0, lambda)){
      d[i] <- 0
    } else if (U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}

# Interpolacao
qpois3<-function(p, lambda, mx=1600){
  q  <- 0:mx
  tt <- p
  pv <- ppois(q, lambda)
  v <- approxfun(pv, q, method="linear") #linear e default
  for(i in 1:length(tt)){
    if(tt[i] <= ppois(0, lambda)){
      tt[i] <- 0
    } else if(tt[i] >= 0.999999999999999999){
      tt[i] <- Inf
    }else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.001)
system.time(q1<-qpois(p, 4))
system.time(q2<-qpois2(p, 4)) 
system.time(q3<-qpois3(p, 4)) 
plot(p, q1)
points(p, q2, col=2)
points(p, q3, col=3)

system.time(q1<-qpois(p, 40))
system.time(q2<-qpois2(p, 40)) 
system.time(q3<-qpois3(p, 40)) 

q<-0:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)
qpois3(p, 4)

q<-3:15
p<-ppois(q, 4)
qpois(p, 4)
qpois2(p, 4)
qpois3(p, 4)

q<-0:150
p<-ppois(q, 40)
qpois(p, 40)
qpois2(p, 40)
qpois3(p, 40)
qpois3(p, 40, mx=150)