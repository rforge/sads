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
library(untb)
setClass("fitsad",representation("mle2", sad="character", trunc="numeric"))

##################################################################################
### Testes da fitlogser
#################################################################################
#Funcao de ajuste para distribuicao log-serie com truncagem
fitlogser2 <- function(x, rich, size, trunc = 0, start.value, ...){
  dots <- list(...)
  if (min(x)<=trunc){
    stop("truncation point should be lower than the lowest data value")
  }
  if (!missing(x)){
    S <- length(x)
    N <- sum(x)
    if (!missing(size)| !missing(rich)){
      warning(paste("Model fitted with size = ", N, " and rich = ", S, " \n calculated from supplied abundances"))
    }
  }
  if (missing(x) & !missing(size) & !missing(rich)){
    S <- rich
    N <- size
  }
  if (missing(x) & missing(size) & missing(rich)){
    stop("Please provide size and species number or a vector of species abundances")
  }
  if (missing(start.value)){
    f1 <- function(a) {
      S + a*log((a/(a + N)))
    }
    sol <- uniroot(f1, interval = c(1/N, N))
    alfa <- sol$root
    X <- N/(N + alfa)
  } else{
    alfa <- start.value
  }
  if (!missing(x)){
    if (trunc >=1){
      LL <- function(alpha) -sum(trunc("dls", x, N, alpha, trunc = trunc, log = T))
    } else{
      LL <- function(alpha) -sum(dls(x, N, alpha, log = T))
    }
    result <- mle2(LL, start = list(alpha = alfa), data = list(x = x), ...)
    new("fitsad", result, sad = "ls", trunc = trunc)
  }
  else new("fitsad", coef = c(alpha = alfa), fullcoef = c(alpha = alfa), sad = "ls", trunc = trunc)
}

set.seed(2012)
## Uma amostra Poisson de uma lognormal
samp1 <- rsad(100, frac = 0.15, sad = lnorm, samp = "Poisson", meanlog = 3, sdlog = 2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N = sum(samp1), S = length(samp1), nmax = sum(samp1))

#Testes
##sem truncagem
x1 <- samp1

nvl1 = function(alpha) -sum(dls(x1, N = sum(x1), alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 1)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2

nvl2 = function(alpha) -sum(dls(x2, N = sum(x2), alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 1)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
####Com chute inicial
(model1 <- fitlogser2(x1, star.value = 1))
attributes(model1)
(model2 <- fitlogser2(x2, star.value = 1))
attributes(model2)
####Sem chute inicial
(model1 <- fitlogser2(x1))
attributes(model1)
(model2 <- fitlogser2(x2))
attributes(model2)
#undebug(fitlogser2)

##com truncagem
x1 <- samp1[samp1 != 1]

nvl1 = function(alpha) -sum(trunc("dls", x1, N = sum(x1), trunc=1, alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 13)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(trunc("dls", x2, N = sum(x2), trunc=1, alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
(tmodel1 <- fitlogser2(x1, trunc = 1, star.value = 13))
attributes(tmodel1)
(tmodel2 <- fitlogser2(x2, trunc = 1, star.value = 13))
attributes(tmodel2)
#debug(fitlogser2)

###############################################################################
####  Graficos e ajustes 
###############################################################################
# sem truncagem
## Amostra 1
x1 <- samp1

nvl1 = function(alpha) -sum(dls(x1, N = sum(x1), alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 1)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

##Amostra 2
x2 <- samp2

nvl2 = function(alpha) -sum(dls(x2, N = sum(x2), alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 1)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
Nm = seq(sum(x1)-200, sum(x1) + 200, 0.1)
alfa = seq(12, 13.5, 0.1)
lls = function(alfa, Nm, x1) -sum(dls(x1, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x1)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=80)
abline(v=12.437, col="red", lty=2)                     
abline(h=sum(x1), col="red", lty=2)                   

##Grafico 2
Nm = seq(sum(x2)-200, sum(x2)+200, 0.1)
alfa = seq(12, 14, 0.1)
lls = function(alfa, Nm, x2) -sum(dls(x2, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x2)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=70)
abline(v=13.1107, col="red", lty=2)                     
abline(h=sum(x2), col="red", lty=2)


# com truncagem em um
##Amostra 1
x1 <- samp1[samp1 != 1]

nvl1 = function(alpha) -sum(trunc("dls", x1, N = sum(x1), trunc=1, alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 15)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)                             

##Amostra 2
x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(trunc("dls", x2, N = sum(x2), trunc=1, alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
Nm = seq(sum(x1)-200, sum(x1)+200, 0.1)
alfa = seq(12, 13, 0.1)
lls = function(alfa, Nm, x1) -sum(trunc("dls", x1, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x1)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=30)
abline(v=12.2329, col="red", lty=2)                     
abline(h=sum(x1), col="red", lty=2)

#Grafico 2
Nm = seq(sum(x2)-200, sum(x2)+200, 0.1)
alfa = seq(11, 14, 0.1)
lls = function(alfa, Nm, x2) -sum(trunc("dls", x2, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x2)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=40)
abline(v=12.1971, col="red", lty=2)                     
abline(h=sum(x2), col="red", lty=2)

#############################################################################
#Metodos otimizacao (poisson lognormal)
############################################################################
fitpoilog2 <- function(x, trunc = 0, ...){
  dots <- list(...)
  if (min(x)<=trunc){
    stop("truncation point should be lower than the lowest data value")
  }
  if (trunc >= 1){
    pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
    LL <- function(mu, sig) -sum(trunc("dpoilog", x, mu, sig, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), ...)
  } else{
    pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
    LL <- function(mu, sig) -sum(trunc("dpoilog", x, mu, sig, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = as.list(pl.par), data = list(x = x), ...)
  }  
  new("fitsad", result, sad="poilog", trunc = trunc)
}

#Testes
##truncada em zero
x1 <- samp1               

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="L-BFGS-B")
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="L-BFGS-B")
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpoilog2
(model1 <- fitpoilog2(x1))
attributes(model1)
(model2 <- fitpoilog2(x2))
attributes(model2)
#undebug(fitpoilog2)

##truncada em um
x1 <- samp1[samp1 != 1]

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="CG")
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="BFGS")
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="CG")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="BFGS")
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
(tmodel1 <- fitpoilog2(x1, trunc = 1, method="CG"))
summary(tmodel1)
(tmodel1 <- fitpoilog2(x1, trunc = 1, method="CG", control=list(maxit = 1000)))
summary(tmodel1)
attributes(tmodel1)
(tmodel2 <- fitpoilog2(x2, trunc = 1))
attributes(tmodel2)
#debug(fitlogser2)


###############################################################################
####  Graficos e ajustes 
###############################################################################
#truncagem em zero
##Amostra 1
x1 <- samp1               

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="SANN")$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="Nelder-Mead") 
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="L-BFGS-B")
summary(nvl1.mle)                            
logLik(nvl1.mle)

##Amostra 2
x2 <- samp2

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="SANN")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="Nelder-Mead")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="L-BFGS-B")
summary(nvl2.mle)
logLik(nvl2.mle)

##Grafico 1
m = seq(0, 5, 0.01)
s = seq(1, 4, 0.1)
lpoilog = function(m, s, x1) -sum(trunc("dpoilog", x1, mu=m, sig=s, trunc=0, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x1)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=1.958, col="red", lty=2)                     
abline(h=1.66, col="red", lty=2)

##Grafico 2
m = seq(0, 5, 0.01)
s = seq(1, 4, 0.1)
lpoilog = function(m, s, x2) -sum(trunc("dpoilog", x2, mu=m, sig=s, trunc=0, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x2)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=2.474, col="red", lty=2)                     
abline(h=1.6778, col="red", lty=2)                   


#truncagem em um
##Amostra 1
x1 <- samp1[samp1 != 1]

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="SANN", control=list(maxit=10000))
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="CG")
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="BFGS")
summary(nvl1.mle)                            
logLik(nvl1.mle)   

##Amostra 2
x2 <- samp2[samp2 != 1]

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="SANN", control=list(maxit=10000)) 
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="CG")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="BFGS")
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
m = seq(-1, 1, 0.01)
s = seq(1, 4, 0.1)
lpoilog = function(m, s, x1) -sum(trunc("dpoilog", x1, trunc=1, mu=m, sig=s, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x1)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=-0.09, col="red", lty=2)                     
abline(h=2.670, col="red", lty=2)                   

##Grafico 2
m = seq(1, 5, 0.01)
s = seq(1, 3, 0.1)
lpoilog = function(m, s, x2) -sum(trunc("dpoilog", x2, trunc=1, mu=m, sig=s, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x2)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=3.00814, col="red", lty=2)                     
abline(h=1.37716, col="red", lty=2)                   

#######################
## Teste 2
#######################
set.seed(42)
## Uma amostra Poisson de uma lognormal
samp1 <- rsad(100, frac = 0.15, sad = lnorm, samp = "Poisson", meanlog = 3, sdlog = 2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N = sum(samp1), S = length(samp1), nmax = sum(samp1))

#Testes
##sem truncagem
x1 <- samp1

nvl1 = function(alpha) -sum(dls(x1, N = sum(x1), alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 1)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2

nvl2 = function(alpha) -sum(dls(x2, N = sum(x2), alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 1)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
####Com chute inicial
(model1 <- fitlogser2(x1, star.value = 1))
attributes(model1)
(model2 <- fitlogser2(x2, star.value = 1))
attributes(model2)
####Sem chute inicial
(model1 <- fitlogser2(x1))
attributes(model1)
(model2 <- fitlogser2(x2))
attributes(model2)
#undebug(fitlogser2)

##com truncagem
x1 <- samp1[samp1 != 1]

nvl1 = function(alpha) -sum(trunc("dls", x1, N = sum(x1), trunc=1, alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 13)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(trunc("dls", x2, N = sum(x2), trunc=1, alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
(tmodel1 <- fitlogser2(x1, trunc = 1, star.value = 13))
attributes(tmodel1)
(tmodel2 <- fitlogser2(x2, trunc = 1, star.value = 13))
attributes(tmodel2)
#undebug(fitlogser2)


###############################################################################
####  Graficos e ajustes 
###############################################################################
# sem truncagem
## Amostra 1
x1 <- samp1

nvl1 = function(alpha) -sum(dls(x1, N = sum(x1), alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 1)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

##Amostra 2
x2 <- samp2

nvl2 = function(alpha) -sum(dls(x2, N = sum(x2), alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 1)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
Nm = seq(sum(x1)-200, sum(x1) + 200, 0.1)
alfa = seq(17, 18, 0.1)
lls = function(alfa, Nm, x1) -sum(dls(x1, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x1)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=80)
abline(v=17.5083, col="red", lty=2)                     
abline(h=sum(x1), col="red", lty=2)                   

##Grafico 2
Nm = seq(sum(x2)-200, sum(x2)+200, 0.1)
alfa = seq(18, 19, 0.1)
lls = function(alfa, Nm, x2) -sum(dls(x2, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x2)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=70)
abline(v=18.6839, col="red", lty=2)                     
abline(h=sum(x2), col="red", lty=2)


# com truncagem em um
##Amostra 1
x1 <- samp1[samp1 != 1]

nvl1 = function(alpha) -sum(trunc("dls", x1, N = sum(x1), trunc=1, alpha, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 15)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)                             

##Amostra 2
x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(trunc("dls", x2, N = sum(x2), trunc=1, alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
Nm = seq(sum(x1)-200, sum(x1)+200, 0.1)
alfa = seq(18, 19, 0.1)
lls = function(alfa, Nm, x1) -sum(trunc("dls", x1, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x1)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=30)
abline(v=18.2431, col="red", lty=2)                     
abline(h=sum(x1), col="red", lty=2)

#Grafico 2
Nm = seq(sum(x2)-200, sum(x2)+200, 0.1)
alfa = seq(18, 19, 0.1)
lls = function(alfa, Nm, x2) -sum(trunc("dls", x2, N = Nm, alpha=alfa, log=T))
llikls = Vectorize(lls, c("alfa", "Nm"))                                   
sup.mat = outer(alfa, Nm, llikls, x2)                                      
contour(alfa, Nm, sup.mat, xlab="alpha", ylab="N", nlevels=40)
abline(v=18.6594, col="red", lty=2)                     
abline(h=sum(x2), col="red", lty=2)

#############################################################################
#Metodos otimizacao (poisson lognormal)
############################################################################
#Testes
##truncada em zero
x1 <- samp1               

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="L-BFGS-B")
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="L-BFGS-B")
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpoilog2
(model1 <- fitpoilog2(x1))
attributes(model1)
(model2 <- fitpoilog2(x2))
attributes(model2)
#undebug(fitpoilog2)

##truncada em um
x1 <- samp1[samp1 != 1]

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="CG")
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="BFGS")
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="CG")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="BFGS")
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
(tmodel1 <- fitpoilog2(x1, trunc = 1))
summary(tmodel1)
attributes(tmodel1)
(tmodel2 <- fitpoilog2(x2, trunc = 1, method="Nelder-Mead"))
attributes(tmodel2)
#debug(fitlogser2)


###############################################################################
####  Graficos e ajustes 
###############################################################################
#truncagem em zero
##Amostra 1
x1 <- samp1               

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="SANN")$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="Nelder-Mead") 
nvl1.mle <- mle2(nvl1, start=as.list(pl.par), control=list(trace=T, maxit=5000), method="L-BFGS-B")
summary(nvl1.mle)                            
logLik(nvl1.mle)

##Amostra 2
x2 <- samp2

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, mu, sig, trunc=0, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="SANN")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="Nelder-Mead")
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="L-BFGS-B")
summary(nvl2.mle)
logLik(nvl2.mle)

##Grafico 1
m = seq(0, 2, 0.01)
s = seq(1, 3, 0.1)
lpoilog = function(m, s, x1) -sum(trunc("dpoilog", x1, mu=m, sig=s, trunc=0, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x1)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=1.20130, col="red", lty=2)                     
abline(h=1.87336, col="red", lty=2)

##Grafico 2
m = seq(0, 2, 0.01)
s = seq(1, 3, 0.1)
lpoilog = function(m, s, x2) -sum(trunc("dpoilog", x2, mu=m, sig=s, trunc=0, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x2)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=1.11004, col="red", lty=2)                     
abline(h=1.80397, col="red", lty=2)                   


#truncagem em um
##Amostra 1
x1 <- samp1[samp1 != 1]

nvl1 = function(mu, sig) -sum(trunc("dpoilog", x1, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="SANN", control=list(maxit=10000))
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="CG")
nvl1.mle = mle2(nvl1, start=as.list(pl.par), method="BFGS")
summary(nvl1.mle)                            
logLik(nvl1.mle)   

##Amostra 2
x2 <- samp2[samp2 != 1]

nvl2 = function(mu, sig) -sum(trunc("dpoilog", x2, trunc=1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="SANN", control=list(maxit=10000)) 
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="CG")
summary(nvl2.mle)                            
logLik(nvl2.mle)

##Grafico 1
m = seq(0, 2, 0.01)
s = seq(1, 3, 0.1)
lpoilog = function(m, s, x1) -sum(trunc("dpoilog", x1, trunc=1, mu=m, sig=s, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x1)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v = 1.49133, col="red", lty=2)                     
abline(h = 1.73427, col="red", lty=2)                   

##Grafico 2
m = seq(1, 2, 0.01)
s = seq(1, 3, 0.1)
lpoilog = function(m, s, x2) -sum(trunc("dpoilog", x2, trunc=1, mu=m, sig=s, log=T))
llikpoilog = Vectorize(lpoilog, c("m","s"))                                   
sup.mat = outer(m, s, llikpoilog, x2)                                      
contour(m, s, sup.mat, xlab="mu", ylab="sig", nlevels=60)
abline(v=1.76150, col="red", lty=2)                     
abline(h=1.47503, col="red", lty=2)             


####################################################################################
####### Ajuste Zipf
####################################################################################
fitzipf <- function(x, trunc = 0, start.value, ...){
  dots <- list(...)
  N <- max(x)
  if (min(x)<=trunc){
    stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    LL <- function(s) -sum(dzipf(x, N, s, log = TRUE))
    sss <- mean(x)
    sss <- mle2(LL, start = as.list(sss), data = list(x = x), ...)
  } else{
    sss <- start.value   #Incluir o chute inicial automatico
  }
  if (trunc >= 1){
    LL <- function(s) -sum(trunc("dzipf", x, N, s, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = as.list(sss), data = list(x = x), ...)
  } else{
    LL <- function(s) -sum(dzipf(x, N, s, log = TRUE))
    result <-  mle2(LL, start = as.list(sss), data = list(x = x), ...)
  }  
  new("fitsad", result, sad="zipf", trunc = trunc)
}

fr <- c(26486, 12053, 5052, 3033, 2536, 2391, 1444, 1220, 1152, 1039)
p <- fr/sum(fr)
lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
opt.f <- function(s) sum((log(p) - lzipf(s, 10))^2)
(opt <- optimize(opt.f, c(0.5, 10)))

LL1 <-function(s) sum(fr*(s*log(1:10)+log(sum(1/(1:10)^s))))
fit <- mle2(LL, start = list(s = 1))
summary(fit)

s.sq <- opt$minimum
s.ll <- fit@coef

plot(1:10,p,log="xy")
lines(1:10, exp(lzipf(s.sq, 10)),col=2)
lines(1:10, exp(lzipf(s.ll, 10)),col=3)

fr <- y
LL <- function(s) -sum(dzipf(fr, N = 10, s, log = TRUE))
result <-  mle2(LL, start = list(s = 1.45), method="Nelder-Mead")
summary(result)

LL2 <- function(s) {
  S <- exp(s)
  -sum(dzipf2(samp1, N = length(samp1), S, log=TRUE))
}
result <- mle2(LL2, start = list(s = log(2))) 
summary(result)   

LL3 <- function(s) -sum(dzipf2(samp1, N = length(samp1), s, log=TRUE))
result <- mle2(LL3, start = list(s = 1.018437)) 
summary(result)

res.prf=profile(result)
LL <- function(s) -sum(dzipf3(samp1, N = length(samp1), s, log=TRUE))
result <- mle2(LL, start = list(s = 2), method = "Brent", upper=0.8, lower=0) 
summary(result)

dzipf2 <- function(x, N, s, log = F) {
  y <- -s*log(x)-log(sum(1/(1:N)^s))
  if(log) return(y)
  else return(exp(y))
}

dzipf3 <- function(x, N, s, log = F){
  y <- 1/(sum(1/(1:N)^s)*x^s)
  if (log) return(log(y))
  else return (y)
}

dzipf(1:10, 10, 1.451385, log=T)
dzipf2(1:10, N =10, 1.451385, log=T)
dzipf3(1:10, 10, 1.451385,log=T)

dzipf(fr, 10, 1.451385, log=T)
dzipf2(fr, 10, 1.451385, log=T)
dzipf3(fr, 10, 1.451385, log=T)

zdata = data.frame(y = 1:5, ofreq = c(63, 14, 5, 1, 2))
fit = vglm(y ~ 1, zipf, zdata, trace = TRUE, weight = ofreq, crit = "coef")
fit = vglm(y ~ 1, zipf(link = identity, init = 3.4), zdata,
           trace = TRUE, weight = ofreq)
fit@misc$N
(shat = Coef(fit))
with(zdata, weighted.mean(y, ofreq))
fitted(fit, matrix = FALSE)
plot(zdata, dzipf(zdata, 5, 2.343904))


