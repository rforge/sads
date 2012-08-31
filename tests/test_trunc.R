library(sads); library(VGAM) #distrb. zipf (dzipf, pzipf)
source("../sads/pkg/sads/R/trunc.R")
source("../sads/pkg/sads/R/dtrunc.R")
source("../sads/pkg/sads/R/ptrunc.R")
source("../sads/pkg/sads/R/qtrunc.R")
source("../sads/pkg/sads/R/dls.R")
source("../sads/pkg/sads/R/pls.R")
source("../sads/pkg/sads/R/ppoilog.R")
source("../sads/pkg/sads/R/qls.R")
source("../sads/pkg/sads/R/qpoilog.R")
source("../sads/pkg/sads/R/qzipf.R")

source("../sads/pkg/sads/R/dpoilog.R")
source("../sads/pkg/sads/R/dpoilog2.R")

source("../sads/pkg/sads/R/ppoilog2.R")
source("../sads/pkg/sads/R/ppoilog4.R")

source("../sads/pkg/sads/R/qpoilog2.R")

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
# Comparar as funcoes de ppoilog
# Troquei os nomes das funcoes ppoilog e ppoilog4(com truncagem em zero)
x<- 0:20
ppoilog(x, mu = 1, sig = 1)
ppoilog2(x, mu = 1, sig = 1)
ppoilog4(x, mu = 1, sig = 1)

trunc("ppoilog", x, mu = 1, sig = 1)
trunc("ppoilog2", x, mu = 1, sig = 1)
trunc("ppoilog4", x, mu = 1, sig = 1)

trunc("ppoilog", x, mu = 1, sig = 1, trunc=0)
trunc("ppoilog2", x, mu = 1, sig = 1, trunc=0)
trunc("ppoilog4", x, mu = 1, sig = 1, trunc=0)

x<-0:100
system.time(trunc("ppoilog2", x, mu = 10, sig = 1, trunc=0))
system.time(trunc("ppoilog4", x, mu = 10, sig = 1, trunc=0))

#Troquei os nomes das funcoes dpoilog e dpoilog2(com truncagem em zero)
x<- 0:20
dpoilog(x, mu = 1, sig = 1)
dtrunc("poilog", x, mu = 1, sig = 1)
trunc("dpoilog", x, mu = 1, sig = 1)

dtrunc("poilog", x, mu = 1, sig = 1, trunc = 0)
trunc("dpoilog", x, mu = 1, sig = 1, trunc = 0)

q<-0:20
ppoilog(q, mu = 1, sig = 1)
ptrunc("poilog", q, mu = 1, sig = 1)
trunc("ppoilog", q, mu = 1, sig = 1)

ptrunc("poilog", q, mu = 1, sig = 1, trunc = 0)
trunc("ppoilog", q, mu = 1, sig = 1, trunc = 0)
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
q<-0:1000
system.time(d1<-ppoilog(q, 10, 13))  # media de 0.012
system.time(d4<-ppoilog4(q, 10, 13)) # media de 5.404

q<-0:100
system.time(d1<-ppoilog(q, 10, 13))  # media de 0.004
system.time(d4<-ppoilog4(q, 10, 13)) # media de 0.064
plot(q, d1)
points(q, d4, col=2, pch="*")

q<-0:10
trunc("ppoilog", q, mu=10, sig=13, trunc=2)
trunc("ppoilog4", q, mu=10, sig=13, trunc=2)


q<-1:1000
system.time(d1<-pls(q, 1000, 1.3))  # Media de 0.000
system.time(d4<-pls4(q, 1000, 1.3)) # Media de 0.184

q<-1:10
system.time(d1<-pls(q, 1000, 1.3))  # media de 0.004
system.time(d4<-pls4(q, 1000, 1.3)) # media de 0.064
plot(q, d1)
points(q, d4, col = 2, pch="*")

q<-1:10
trunc("pls", q, N = 1000, alpha = 1.3, trunc = 2)
trunc("pls4", q, N = 1000, alpha = 1.3, trunc = 2)

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

##############################################################
# Testes quantilicos para log-serie e poisson-lognormal     ##
##############################################################
#Log-serie
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

qls2<-function(p, N, alpha, mx=1600){
  q  <- 1:mx
  tt <- p
  pv <- pls(q, N, alpha)
  v <- approxfun(pv, q, method="linear") #linear e default
  for(i in 1:length(tt)){
    if(tt[i] <= pls(1, N, alpha)){
      tt[i] <- 1
    } else if(tt[i] >= 0.999999999999999999){
      tt[i] <- Inf
    }else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.001)
system.time(q1<-qls(p, 10, 2)) #media 0.480
system.time(q2<-qls2(p, 10, 2)) #media 0.060
plot(p, q1)
points(p, q2, col=2)

system.time(q1<-qls(p, 100, 2)) #media 6.601
system.time(q2<-qls2(p, 100, 2)) #media 0.112
plot(p, q1)
points(p, q2, col=2)

q<-1:15
p<-pls(q, 10, 2)
qls(p, 10, 2)
qls2(p, 10, 2)

q<-3:15
p<-pls(q, 10, 2)
qls(p, 10, 2)
qls2(p, 10, 2)

q<-1:150
p<-pls(q, 100, 2)
qls(p, 100, 2)
qls2(p, 100, 2)
qls2(p, 100, 2, mx=100)

# Poisson-lognormal
qpoilog <- function(p, mu, sig){
  d <- NULL
  busca <- function(U1, U2){
    repeat{
      tt <- ppoilog(U2, mu, sig)
      U2 <- ifelse(tt >= U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt < U1, U2+1, U2)
      a2 <- U2
      if (ppoilog(min(a1, a2), mu, sig) < U1 & U1 <= ppoilog(max(a1, a2), mu, sig)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1000, max=3000000))
    if(U1 <= ppoilog(0, mu, sig)){
      d[i] <- 0
    } else if(U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}

qpoilog2 <- function(p, mu, sig, mx = 115000){
  q  <- 0:mx
  tt <- p
  pv <- ppoilog(q, mu, sig)
  v <- approxfun(pv, q, method="linear") #linear e default
  for(i in 1:length(tt)){
    if (tt[i] <= ppoilog(0, mu, sig)){
      tt[i] <- 0
    } else if(tt[i] >= 0.999999999999999999){
      tt[i] <- Inf
    } else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.01)
system.time(q1<-qpoilog(p, 1, 2)) #media 77.113
system.time(q2<-qpoilog2(p, 1, 2)) #media 0.172
plot(p, q1)
points(p, q2, col=2, pch="*")

system.time(q1<-qpoilog(p, 10, 2)) #media 6.601
system.time(q2<-qpoilog2(p, 10, 2)) #media 0.112
plot(p, q2)
points(p, q2, col=2)

q<-0:15
p<-ppoilog(q, 10, 2)
qpoilog(p, 10, 2)
qpoilog2(p, 10, 2)

q<-3:15
p<-ppoilog(q, 10, 2)
qpoilog(p, 10, 2)
qpoilog2(p, 10, 2)

q<-1:150
p<-ppoilog(q, 100, 2)
qpoilog(p, 100, 2)
qpoilog2(p, 100, 2)
qpoilog2(p, 100, 2, mx=100)

#Power
qpower2<-function(p, s){
  d<-NULL
  busca <- function(U1, U2){
    repeat{
      tt <- ppower(U2, s)
      U2 <- ifelse(tt>=U1, U2-1, U2)
      a1 <- U2
      U2 <- ifelse(tt<U1, U2+1, U2)
      a2 <- U2
      if(ppower(min(a1, a2), s) < U1 & U1 <= ppower(max(a1, a2), s)){
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)){
    U1 <- p[i]
    U2 <- round(runif(1, min=1, max=length(p)))
    if(U1 <= ppower(1, s)){
      d[i] <- 1
    } else if (U1 >= 0.999999999999999999){
      d[i] <- Inf
    } else{
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}

qpower3<-function(p, s, mx=1600){
  q  <- 1:mx
  tt <- p
  pv <- ppower(q, s)
  v <- approxfun(pv, q, method="linear") #linear e default
  for(i in 1:length(tt)){
    if(tt[i] <= ppower(1, s)){
      tt[i] <- 1
    } else if(tt[i] >= 0.999999999999999999){
      tt[i] <- Inf
    }else {
      tt[i] <- ceiling(v(tt[i]))
    }
  }
  return(tt)
}

p<-seq(0, 1, 0.01)
system.time(q0<-qpower(p, 2)) #media 2.548
system.time(q1<-qpower2(p, 2)) #media 2.756
system.time(q2<-qpower3(p, 2)) #media 0.624
plot(p, q0)
points(p, q1, col=2)
points(p, q2, col=3)

system.time(q1<-qls(p, 100, 2)) #media 6.601
system.time(q2<-qls2(p, 100, 2)) #media 0.112
plot(p, q1)
points(p, q2, col=2)

q<-1:15
p<-ppower(q, 2)
qpower(p, 2)
qpower2(p, 2)
qpower3(p, 2)

q<-3:15
p<-ppower(q, 2)
qpower(p, 2)
qpower2(p, 2)
qpower3(p, 2)

q<-1:150
p<-ppower(q, 2)
qpower(p, 2)
qpower2(p, 2)
qpower3(p, 2)


###############################################################################
## Testes Extras
###############################################################################
#p-functions
q <- 1:20
pls(q, 100, 3.5)
pls(q, 100, 3.5, log = T)
log(pls(q, 100, 3.5))

ppoilog(q, 1, 2)
ppoilog(q, 1, 2, log = T)
log(ppoilog(q, 1, 2))

ppower(q, 2)
ppower(q, 2, log = T)
log(ppower(q, 2))

pzipf(q, 20, 2)
pzipf(q, 20, 2, log = T)
log(pzipf(q, 20, 2))

#q-functions
q <- 1:20
p0 <- pls(q, 100, 3.5)
qls(p0, 100, 3.5)
p1 <- pls(q, 100, 3.5, log = T)
qls(p1, 100, 3.5, log = T)
p0 <- pls(q, 100, 3.5, lower = F)
qls(p0, 100, 3.5, lower = F)
p1 <- pls(q, 100, 3.5, lower = F, log = T)
qls(p1, 100, 3.5, lower = F, log = T)

p0 <- ppoilog(q, 1, 2)
qpoilog(p0, 1, 2)
p1 <- ppoilog(q, 1, 2, log = T)
qpoilog(p1, 1, 2, log = T)
p0 <- ppoilog(q, 1, 2, lower = F)
qpoilog(p0, 1, 2, lower = F)
p1 <- ppoilog(q, 1, 2, lower = F, log = T)
qpoilog(p1, 1, 2, lower = F, log = T)

p0 <- ppower(q, 2)
qpower(p0, 2)
p1 <- ppower(q, 2, log = T)
qpower(p1, 2, log = T)
p0 <- ppower(q, 2, lower = F)
qpower(p0, 2, lower = F)
p1 <- ppower(q, 2, lower = F, log = T)
qpower(p1, 2, lower = F, log = T)

p0 <- pzipf(q, 20, 2)
qzipf(p0, 20, 2)
p1 <- pzipf(q, 20, 2, log = T)
qzipf(p1, 20, 2, log = T)
p0 <- pzipf(q, 20, 2, lower = F)
qzipf(p0, 20, 2, lower = F)
p1 <- pzipf(q, 20, 2, lower = F, log = T)
qzipf(p1, 20, 2, lower = F, log = T)


x <- 0:10
dtrunc("pois", x, coef = 2)
dpois(x, 2) 
p1 <- dtrunc("pois", x, coef = 2, log=T)
dpois(x, 2, log = T)

(q0 <- ptrunc("pois", x, coef = 2))
ppois(x, 2)
(q1 <- ptrunc("pois", x, coef = 2, log = T))
ppois(x, 2, log = T)

qtrunc("pois", q0, lambda = 2)
qtrunc("pois", q1, lambda = 2)
qtrunc("pois", q0, lambda = 2, log = T)
qtrunc("pois", q1, lambda = 2, log = T)


#################################################
### Truncagem
#################################################
x <- 1:15
dtrunc("pois", x, coef = 2)
dpois(x, 2)
dtrunc("pois", x, coef = 2, log = T)
dpois(x, 2, log = T)
ptrunc("pois", x, coef = 2)
ppois(x, 2)
ptrunc("pois", x, coef = 2, log = T)
ppois(x, 2, log = T)
dtrunc("pois", x, coef = 2, trunc = 2)
dtrunc2("pois", x, lambda = 2, trunc = 2)
dtrunc("pois", x, coef = 2, trunc = 2, log = T)
dtrunc2("pois", x, lambda = 2, trunc = 2, log = T)
ptrunc("pois", x, coef = 2, trunc = 2)
ptrunc2("pois", x, lambda = 2, trunc = 2)
ptrunc("pois", x, coef = 2, trunc = 2, log = T)
ptrunc2("pois", x, lambda = 2, trunc = 2, log = T)

x <- seq(-4, 4, 0.5)
dtrunc("norm", x, coef = list(2, 3))
dtrunc2("norm", x, 2, 3)
dnorm(x, 2, 3)
dtrunc("norm", x, coef = list(2, 3), log = T)
dnorm(x, 2, 3, log = T)
ptrunc("norm", x, coef = list(2, 3))
pnorm(x, 2, 3)
ptrunc("norm", x, coef = list(2, 3), log = T)
pnorm(x, 2, 3, log = T)
dtrunc("norm", x, coef = list(2, 3), trunc = -2)
dtrunc2("norm", x, 2, 3, trunc = -2)
dtrunc("norm", x, coef = list(2, 3), trunc = -2, log = T)
dtrunc2("norm", x, 2, 3, trunc = -2, log = T)
ptrunc("norm", x, coef = list(2, 3), trunc = -2)
ptrunc2("norm", x, 2, 3, trunc = -2)
ptrunc("norm", x, coef = list(2, 3), trunc = -2, log = T)
ptrunc2("norm", x, 2, 3, trunc = -2, log = T)


x <- 1:15
p1 <- ptrunc("pois", x, coef = 2, lower = T, log = F)
qtrunc("pois", p1, coef = 2, lower = T, log = F)
p2 <- ptrunc("pois", x, coef = 2, lower = F, log = F)
qtrunc("pois", p2, coef = 2, lower = F, log = F)
p3 <- ptrunc("pois", x, coef = 2, lower = T, log = T)
qtrunc("pois", p3, coef = 2, lower = T, log = T)
p4 <- ptrunc("pois", x, coef = 2, lower = F, log = T)
qtrunc("pois", p4, coef = 2, lower = F, log = T)

x <- 1:15
p1 <- ptrunc("pois", x, coef = 2, trunc = 3, lower = T, log = F)
qtrunc("pois", p1, coef = 2, trunc = 3, lower = T, log = F)
p2 <- ptrunc("pois", x, coef = 2, trunc = 3, lower = F, log = F) #prob
qtrunc("pois", p2, coef = 2, trunc = 3, lower = F, log = F)
p3 <- ptrunc("pois", x, coef = 2, trunc = 3, lower = T, log = T)
qtrunc("pois", p3, coef = 2, trunc = 3, lower = T, log = T)
p4 <- ptrunc("pois", x, coef = 2, trunc = 3, lower = F, log = T) #prob
qtrunc("pois", p4, coef = 2, trunc = 3, lower = F, log = T)

x <- seq(-4, 4, 0.5)
p1 <- ptrunc("norm", x, coef = list(2, 3), trunc = -3, lower = T, log = F)
qtrunc("norm", p1, coef = list(2, 3), trunc = -3, lower = T, log = F)
p2 <- ptrunc("norm", x, coef = list(2, 3), trunc = -3, lower = F, log = F) #prob
qtrunc("norm", p2, coef = list(2, 3), trunc = -3, lower = F, log = F)
p3 <- ptrunc("norm", x, coef = list(2, 3), trunc = -3, lower = T, log = T)
qtrunc("norm", p3, coef = list(2, 3), trunc = -3, lower = T, log = T)
p4 <- ptrunc("norm", x, coef = list(2, 3), trunc = -3, lower = F, log = T) #prob
qtrunc("norm", p4, coef = list(2, 3), trunc = -3, lower = F, log = T)