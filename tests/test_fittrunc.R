library(sads)
library(VGAM) #distrb. zipf (dzipf, pzipf)
library(untb)
source("../sads/pkg/sads/R/trunc.R")
source("../sads/pkg/sads/R/dtrunc.R")
source("../sads/pkg/sads/R/ptrunc.R")
source("../sads/pkg/sads/R/qtrunc.R")
source("../sads/pkg/sads/R/dls.R")
source("../sads/pkg/sads/R/pls.R")
source("../sads/pkg/sads/R/qls.R")
source("../sads/pkg/sads/R/fitls.R")
source("../sads/pkg/sads/R/dpoilog.R")
source("../sads/pkg/sads/R/ppoilog.R")
source("../sads/pkg/sads/R/qpoilog.R")
source("../sads/pkg/sads/R/fitpoilog.R")
source("../sads/pkg/sads/R/dzipf.R")
source("../sads/pkg/sads/R/pzipf.R")
source("../sads/pkg/sads/R/qzipf.R")
source("../sads/pkg/sads/R/fitzipf.R")
source("../sads/pkg/sads/R/dpower.R")
source("../sads/pkg/sads/R/ppower.R")
source("../sads/pkg/sads/R/qpower.R")
source("../sads/pkg/sads/R/fitpower.R")

#source("../sads/pkg/sads/R/dpoilog2.R")
#source("../sads/pkg/sads/R/ppoilog2.R")
#source("../sads/pkg/sads/R/ppoilog4.R")
#source("../sads/pkg/sads/R/qpoilog2.R")
setClass("fitsad", representation("mle2", sad="character", trunc="numeric"))

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
(model1 <- fitlogser(x1, start.value = 1))
attributes(model1)
(model2 <- fitlogser(x2, start.value = 1))
attributes(model2)
####Sem chute inicial
(model1 <- fitlogser(x1))
attributes(model1)
(model2 <- fitlogser(x2))
attributes(model2)
#undebug(fitlogser2)

##com truncagem
x1 <- samp1[samp1 != 1]

nvl1 = function(alpha) -sum(dtrunc("ls", x1, coef = as.list(N = sum(x1), alpha = alpha), trunc=1, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 13), data = list(x = x)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(trunc("dls", x2, N = sum(x2), trunc=1, alpha, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitlogser2
(tmodel1 <- fitlogser(x1, trunc = 1, start.value = 13))
attributes(tmodel1)
(tmodel2 <- fitlogser(x2, trunc = 1, start.value = 13))
attributes(tmodel2)
#debug(fitlogser)

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
(model1 <- fitlogser2(x1, start.value = 1))
attributes(model1)
(model2 <- fitlogser2(x2, start.value = 1))
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
(model1 <- fitpoilog(x1))
attributes(model1)
(model2 <- fitpoilog(x2))
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
(tmodel1 <- fitpoilog(x1, trunc = 1))
summary(tmodel1)
attributes(tmodel1)
(tmodel2 <- fitpoilog(x2, trunc = 1, method="Nelder-Mead"))
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
dzipf <- function(x, N, s, log = F) {
  if (any(x < 1)) warning("the zipf's distribution is not set to zero")
  if (s <= 0) stop("s must be greater than zero")
  if (N < 1) stop("N must be positive integer")
  if (!any(is.wholenumber(x))) warning("x must be integer")
  y <- NULL
  for (i in 1:length(x)){
    if(!is.wholenumber(x[i])) y[i] <- -Inf
    else y[i] <- -s*log(x[i])-log(sum(1/(1:N)^s))
  }
  if(log) return(y)
  else return(exp(y))
}

pzipf <- function(q, N, s){
  if (s <= 0) stop("s must be greater than zero")
  if (N < 1) stop("N must be positive integer")
  y <- NULL
  for (i in 1:length(q)){
    y[i] <- log(sum(1/(1:q[i])^s)) - log(sum(1/(1:N)^s))
  }
  return(exp(y))
}

fitzipf <- function(x, N, trunc = 0, start.value, upper = 20, ...){
  if (min(x)<=trunc){
    stop("truncation point should be lower than the lowest data value")
  }
  if(missing(N)){
    N <- length(x)
  }
  if(missing(start.value)){
    p <- x/sum(x)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    print(opt)
    sss <- opt$minimum
    #LL1 <-function(s) sum(x*(s*log(1:N)+log(sum(1/(1:N)^s))))
    #fit <- mle2(LL1, start = list(s = sss))
    #sss <- exp(fit@coef)
    #print(fit)
  }else{
    sss <- start.value
  }
  if(trunc >= 1){
    LL <- function(s) -sum(trunc("dzipf", x, N, s, trunc = trunc, log = TRUE))
    result <-  mle2(LL, start = list(s = sss), data = list(x = x), method = "Brent", lower = 0, upper = upper, ...)
  } else{
    LL <-function(s) sum(x*(((1:N)^s)*(sum(1/(1:N)^s))))
    result <- mle2(LL, start = list(s = sss), data = list(x = x), method="Brent", lower = 0, upper = upper, ...)
    #LL <- function(s) -sum(dzipf(x, N, s, log = TRUE))
    #result <- mle2(LL, start = list(s = sss), data = list(x = x), method="Brent", lower = 0, upper = upper, ...)
  }
  if(abs(as.numeric(result@coef) - upper) < 0.001) warning("Check the upper limit of the function")
  new("fitsad", result, sad="zipf", trunc = trunc)
}

lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
model1 <- fitzipf(samp1)
summary(model2)
p <- samp1/sum(samp1)
plot(1:length(p), sort(p, d=T), log="xy")
lines(1:length(p), exp(lzipf(1.480988, length(p))), col = 2)
lines(1:length(p), exp(lzipf(0.954095, length(p))), col = 3)
plot(profile(model1))

model2 <- fitzipf(x1)
summary(model2)
p <- x1/sum(x1)
plot(1:length(p), sort(p, d = T))
lines(1:length(p), exp(lzipf(1.335487, length(p))), col = 2)
lines(1:length(p), exp(lzipf(0.616874, length(p))), col = 3)
plot(profile(model2))

model3 <- fitzipf(x1, trunc = 1)
summary(model3)
p <- x1/sum(x1)
plot(1:length(p), sort(p, d=T))
lines(1:length(p), exp(lzipf(1.335487, length(p))), col = 2)
lines(1:length(p), exp(lzipf(0.86503, length(p))), col = 3)
plot(profile(model3))


################################## Outros testes
LL2 <- function(s) {
  S <- exp(s)
  -sum(dzipf2(samp1, N = length(samp1), S, log=TRUE))
}
result <- mle2(LL2, start = list(s = log(2))) 
summary(result)   

LL <- function(s) -sum(dzipf3(samp1, N = length(samp1), s, log=TRUE))
result <- mle2(LL, start = list(s = 2), method = "Brent", upper=1.8, lower=0.00) 
summary(result)
res.prf<-profile(result)
plot(res.prf)

p <- samp1/sum(samp1)
lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
(opt <- optimize(opt.f, c(0.5, length(p))))

LL1 <-function(s) sum(samp1*(s*log(1:length(samp1))+log(sum(1/(1:length(samp1))^s))))
fit <- mle2(LL1, start = list(s = 1))
summary(fit)
logLik(fit)

s.sq <- opt$minimum
s.ll <- fit@coef

plot(1:length(p), sort(p, d=T), log="xy")
lines(1:length(p), exp(lzipf(s.sq, length(p))),col=2)
lines(1:length(p), exp(lzipf(s.ll, length(p))),col=3)


fr <- y
LL <- function(s) -sum(dzipf(fr, N = 10, s, log = TRUE))
result <-  mle2(LL, start = list(s = 1.45), method="Nelder-Mead")
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


###############################################################################################
### Log Normal Distribution
###############################################################################################

LL <- function(meanlog, sdlog) - sum(dlnorm(samp1, meanlog, sdlog, log = T))
meanlog <- sum(log(samp1))/length(samp1)
sdlog <- sqrt((sum(log(samp1)-meanlog)^2)/length(samp1))
result <- mle2(LL, start = list(meanlog=exp(meanlog), sdlog= exp(sdlog)), method="SANN")
result <- mle2(LL, start = as.list(result@coef))
summary(result)
result <- mle2(LL, start = list(meanlog=exp(meanlog), sdlog= exp(sdlog)), method="Nel")
summary(result)

fitlnorm <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    meanlog <- sum(log(x))/length(x)
    sdlog <- sqrt((sum(log(x)-meanlog)^2)/length(x))
  } else{
    meanlog <- start.value[1]
    sdlog <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(meanlog, sdlog) -sum(dlnorm(x, meanlog, sdlog, log = TRUE))
  } else {
    LL <- function(meanlog, sdlog) -sum(trunc("dlnorm", x, meanlog, sdlog, trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(meanlog = meanlog, sdlog = sdlog), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="lnorm", trunc = ifelse(missing(trunc), NaN, trunc)) 
}

fitlnorm(samp1)
attributes(fitlnorm(samp1, start=c(1.944075, 1.648148)))
fitlnorm(samp1, trunc = 0.5)
attributes(fitlnorm(samp1, start=c(1.294512, 2.104622), trunc = 0.5))
fitlnorm(samp1, start=c(1.294512, 2.104622), trunc = 0.5, control = list(maxit=1000))
fitlnorm(x1, trunc = 1)
fitlnorm(x1, start=c(1.750650, 1.948772), trunc = 1)
fitlnorm(x1, start=c(1.750650, 1.948772), trunc = 1, control = list(maxit=1000))

###############################################################################################
### Negative Binomial Distribution
###############################################################################################
LL1 <- function(prob) -sum(dnbinom(samp1, size = length(samp1), prob, log = T))
LL2 <- function(mu) -sum(dnbinom(samp1, size = length(samp1), mu, log = T))
size <- length(samp1)
prob <- size/(size+mu)
result1 <- mle2(LL1, start = list(prob = prob), method="Brent", upper = 1, lower = 0)
summary(result1)
mu <- size*(1-as.numeric(result1@coef))/as.numeric(result1@coef)
mu <-mean(samp1)
result2 <- mle2(LL2, start = list(mu = mu), method="CG")
summary(result2)
result <- mle2(LL, start = as.list(result@coef))
summary(result)
optimize(LL2, c(1000, 10000))
result <- mle2(LL, start = list(meanlog=exp(meanlog), sdlog= exp(sdlog)), method="Nel")


fitnbinom <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    LL <- function(prob) -sum(dnbinom(x, size = length(x), prob, log = TRUE))
    phat<- length(x)/(length(x)+mean(x))
    phat <- mle2(LL, start = list(prob = phat), method="Brent", upper = 1, lower = 0)@coef
    mu <- length(x)*(1 - phat)/phat
  } else{
    mu <- start.value
  }
  if (missing(trunc)){
    LL <- function(mu) -sum(dnbinom(x, size = length(x), mu, log = TRUE))
  } else{
    LL <- function(mu) -sum(trunc("dnbinom", x, size = length(x), mu, trunc = trunc, log = TRUE))
  }
  result <- mle2(LL, start = list(mu = mu), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="nbinom", trunc = ifelse(missing(trunc), NaN, trunc))
}

fitnbinom(samp1)
fitnbinom(samp1, start = mean(samp1))
fitnbinom(samp1, trunc = 0)
fitnbinom(samp1, start = mean(samp1), trunc = 0)
fitnbinom(samp1, start = mean(samp1), trunc = 0, control = list(maxit=1000))
fitnbinom(x1, trunc = 1)
fitnbinom(x1, start = mean(x1), trunc = 1)
fitnbinom(x1, start = mean(x1), trunc = 1, control = list(maxit=1000))

##############################################################################################
### Gamma Distribution
###############################################################################################
fitgamma <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    ka <- mean(x)/sd(x)
    theta <- 0.5/(log(mean(x))-mean(log(x)))
  } else{
    ka <- start.value[1]
    theta <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(shape, scale) -sum(dgamma(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(trunc("dgamma", x, shape, scale, trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(shape = ka, scale = theta), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="gamma", trunc = ifelse(missing(trunc), NaN, trunc)) 
}

ka[1] <- (mean(samp1)/sd(samp1))^2
theta[1] <- var(samp1)/mean(samp1)
for(i in 2:100){
  theta[i] <- mean(samp1)/ka[i-1]
  ka[i] <- exp(log(mean(samp1)/(prod(samp1)^(1/length(samp1)))) + digamma(ka[i-1]))
}

thetahat <- (mean(samp1^2) - (mean(samp1))^2)/mean(samp1)
kahat <- (mean(samp1)^2)/(mean(samp1^2) - (mean(samp1)^2))

kahatfunc <- function(ka, xvec){
  n <- length(xvec)
  eq <- -n*digamma(ka)-n*log(mean(xvec))+n*log(ka)+sum(log(xvec));
  eq;
}

kahatfunc(0.06215943, samp1)
kahatfunc(791, samp1)
karoot <- uniroot(kahatfunc, interval=c(.0621,791), xvec=samp1)

ka <- karoot$root
theta <- mean(samp1)/ka

fitgamma(samp1)
attributes(fitgamma(samp1, start = c(0.345477, 0.0070250)))
attributes(fitgamma(samp1, trunc = 0))
fitgamma(samp1, start = c(0.345477, 0.0070250), trunc = 0)
fitgamma(samp1, start = c(0.345477, 0.0070250), trunc = 0, control = list(maxit=1000))
fitgamma(x1, trunc = 1)
fitgamma(x1, start = c(10, 10), trunc = 1)
fitgamma(x1, start = mean(x1), trunc = 1, control = list(maxit=1000))

##############################################################################################
### Weibull Distribution
###############################################################################################
fitweibull <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    ka <- 1
    theta <- mean(x)
    for(i in 1:100){
      theta <- (sum(x^ka)/length(x))^(1/ka)
      ka <- length(x)/(sum(x^ka * log(x)) - sum(log(x))/theta)
    }
  } else{
    ka <- start.value[1]
    theta <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(shape, scale) -sum(dweibull(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(trunc("dweibull", x, shape, scale, trunc = trunc, log = TRUE))
  }  
  result <- mle2(LL, start = list(shape = ka, scale = theta), method="SANN")
  result <- mle2(LL, start = as.list(result@coef), data = list(x = x), ...)
  new("fitsad", result, sad="weibull", trunc = ifelse(missing(trunc), NaN, trunc)) 
}

ka[1] <- 1
theta[1] <- mean(samp1)
for(i in 2:100){
  theta[i] <- (sum(samp1^ka[i-1])/length(samp1))^(1/ka[i-1])
  ka[i] <- length(samp1)/(sum(samp1^ka[i-1] * log(samp1)) - sum(log(samp1))/theta[i-1])
}

fitweibull(samp1)
attributes(fitweibull(samp1, start = c(0.5088313, 17.0280350)))
attributes(fitweibull(samp1, trunc = 0.5))
fitweibull(samp1, start = c(0.5088554, 17.0360020), trunc = 0.5)
fitweibull(samp1, start = c(0.2821154, 1.58950230), trunc = 0.5, control = list(maxit=1000))
fitweibull(x1, trunc = 1)
fitweibull(x1, start = c(0.2893861, 2.0628751), trunc = 1)
fitweibull(x1, start = c(0.2893861, 2.0628751), trunc = 1, control = list(maxit=1000))

fitweibull(samp2)
attributes(fitweibull(samp2, start = c(0.6295015, 26.4126025)))
attributes(fitweibull(samp2, trunc = 0.5))
fitweibull(samp2, start = c(0.4582524, 13.9512424), trunc = 0.5)
fitweibull(samp2, start = c(0.4582421, 13.9512404), trunc = 0.5, control = list(maxit=1000))
fitweibull(x2, trunc = 1)
fitweibull(x2, start = c(0.6955263, 35.6361214), trunc = 1)
fitweibull(x2, start = c(0.6955263, 35.6361214), trunc = 1, control = list(maxit=1000))


##############################################################################################
### Geometric Distribution
###############################################################################################
fitgeom <- function(x, trunc, start.value, ...){
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    phat <- 1/(1 + mean(x))
  } else{
    phat <- start.value
  }
  if (missing(trunc)){
    LL <- function(prob) -sum(dgeom(x, prob, log = T))
  } else{
    LL <- function(prob) -sum(trunc("dgeom", x, prob, trunc = trunc, log = T))
  }
  result <- mle2(LL, start = list(prob = phat), data = list(x = x), method = "Brent", lower = 0, upper = 1, ...)
  new("fitsad", result, sad = "geom", trunc = ifelse(missing(trunc), NaN, trunc))
}

fitgeom(samp1)
attributes(fitgeom(samp1, start = 0.0280350))
attributes(fitgeom(samp1, trunc = 0))
fitgeom(samp1, start = 0.5088554, trunc = 0)
fitgeom(samp1, start = 0.2821154, trunc = 0, control = list(maxit=1000))
fitgeom(x1, trunc = 1)
fitgeom(x1, start = 0.2893861, trunc = 1)
fitgeom(x1, start = 0.2893861, trunc = 1, control = list(maxit=1000))

fitgeom(samp2)
attributes(fitgeom(samp2, start = 0.02516176))
attributes(fitgeom(samp2, trunc = 0))
fitgeom(samp2, start = 0.5088554, trunc = 0)
fitgeom(samp2, start = 0.2821154, trunc = 0, control = list(maxit=1000))
attributes(fitgeom(x2, trunc = 1))
fitgeom(x2, start = 0.2893861, trunc = 1)
fitgeom(x2, start = 0.2893861, trunc = 1, control = list(maxit=1000))

##############################################################################################
### Power-law/Pareto Distribution
###############################################################################################
dpower <- function(x, s, log = FALSE){
  if (any(x < 1)) warning("the zipf's distribution is not set to zero")
  if (s <= 0) stop("s must be greater than zero")
  if (!any(is.wholenumber(x))) warning("x must be integer")
  y <- NULL
  for (i in 1:length(x)){
    if(!is.wholenumber(x[i])) y[i] <- -Inf
    else y[i] <- -s*log(x[i])-log(zeta(s))
  }
  if(log) return(y)
  else return(exp(y))
}

dpower(1:10, 2.34)
dzeta(1:10, 2.34)
dpower(1:10, 3.34)

dzeta(1:10, 2.34, log = T)
dpower(1:10, 3.34, log = T)

ppower <- function(q, s){
  if (s <= 0) stop("s must be greater than zero")
  y <- NULL
  for (i in 1:length(q)){
    y[i] <- log(sum(1/(1:q[i])^s)) - log(zeta(s))
  }
  return(exp(y))
}

plot(1:10, ppower(1:10, 2.34), ylim=c(0, 1))
lines(1:10, dpower(1:10, 2.34), t="l", col=2)

fitpower <- function(x, trunc, start.value, upper = 20, ...){
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    shat <- 2
  } else{
    shat <- start.value
  }
  if (missing(trunc)){
    LL <- function(s) -sum(dpower(x, s, log = T))
  } else{
    LL <- function(s) -sum(trunc("dpower", x, s, trunc = trunc, log = T))
  }
  result <- mle2(LL, start = list(s = shat), data = list(x = x), method = "Brent", lower = 1, upper = upper, ...)
  if(abs(as.numeric(result@coef) - upper) < 0.0000001) warning("check the upper limit of the function")
  new("fitsad", result, sad = "power", trunc = ifelse(missing(trunc), NaN, trunc))
}

fitpower(samp1)
attributes(fitpower(samp1, start = 0.0280350, upper = 1.4))
fitpower(x1, trunc = 1)
attributes(fitpower(x1, start = 1.2893861, trunc = 1))
fitpower(samp1, start = 0.0280350, upper = 1.4)

fitpower(samp2)
attributes(fitpower(samp2, start = 1.0280350))
fitpower(x2, trunc = 1)
attributes(fitpower(x2, start = 1.2893861, trunc = 1))