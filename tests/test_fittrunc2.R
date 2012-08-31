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

set.seed(2012)
## Uma amostra Poisson de uma lognormal
samp1 <- rsad(100, frac = 0.15, sad = lnorm, samp = "Poisson", meanlog = 3, sdlog = 2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N = sum(samp1), S = length(samp1), nmax = sum(samp1))
##################################################################################
### Testes da fitls
#################################################################################
#Funcao de ajuste para distribuicao log-serie com truncagem
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

###teste com funcao fitls
####Com chute inicial
(model1 <- fitls(x1, start.value = 1))
attributes(model1)
(model2 <- fitls(x2, start.value = 1))
attributes(model2)
####Sem chute inicial
(model1 <- fitls(x1))
attributes(model1)
(model2 <- fitls(x2))
attributes(model2)

##com truncagem
x1 <- samp1[samp1 != 1]
nvl1 = function(alpha) -sum(dtrunc("ls", x1, coef = list(N = sum(x1), alpha = alpha), trunc = 1, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 13), data = list(x = x)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(alpha) -sum(dtrunc("ls", x2, coef = list(N = sum(x2), alpha = alpha), trunc=1, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitls
(tmodel1 <- fitls(x1, trunc = 1, start.value = 13))
attributes(tmodel1)
(tmodel2 <- fitls(x2, trunc = 1, start.value = 13))
attributes(tmodel2)

set.seed(42)
## Uma amostra Poisson de uma lognormal
samp1 <- rsad(100, frac = 0.15, sad = lnorm, samp = "Poisson", meanlog = 3, sdlog = 2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N = sum(samp1), S = length(samp1), nmax = sum(samp1))

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

###teste com funcao fitls
####Com chute inicial
(model1 <- fitls(x1, start.value = 1))
attributes(model1)
(model2 <- fitls(x2, start.value = 1))
attributes(model2)
####Sem chute inicial
(model1 <- fitls(x1))
attributes(model1)
(model2 <- fitls(x2))
attributes(model2)

##com truncagem
x1 <- samp1[samp1 != 1]
nvl1 = function(alpha) -sum(dtrunc("ls", x1, coef = list(N = sum(x1), alpha = alpha), trunc = 1, log=TRUE))
nvl1.mle = mle2(nvl1, start=list(alpha = 13), data = list(x = x)) 
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]

nvl2 = function(alpha) -sum(dtrunc("ls", x2, coef = list(N = sum(x2), alpha = alpha), trunc=1, log=TRUE))
nvl2.mle = mle2(nvl2, start=list(alpha = 13)) 
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitls
(tmodel1 <- fitls(x1, trunc = 1, start.value = 13))
attributes(tmodel1)
(tmodel2 <- fitls(x2, trunc = 1, start.value = 13))
attributes(tmodel2)

##################################################################################
### Testes da fitpoilog
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(mu, sig) -sum(dpoilog(x1, mu, sig, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))), zTrunc = F)$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(mu, sig) -sum(dpoilog(x2, mu, sig, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))), zTrunc = F)$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpoilog
(model1 <- fitpoilog(x1))
attributes(model1)
(model2 <- fitpoilog(x2))
attributes(model2)

##truncada em zero
x1 <- samp1               
nvl1 = function(mu, sig) -sum(dtrunc("poilog", x1, coef = list(mu, sig), trunc = 0, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle <- mle2(nvl1, start=as.list(pl.par))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(mu, sig) -sum(dtrunc("poilog", x2, coef = list(mu, sig), trunc = 0, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="L-BFGS-B")
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpoilog
(model1 <- fitpoilog(x1, trunc = 0))
attributes(model1)
(model2 <- fitpoilog(x2, trunc = 0))
attributes(model2)
#undebug(fitpoilog2)

##truncada em um
x1 <- samp1[samp1 != 1]
nvl1 = function(mu, sig) -sum(dtrunc("poilog", x1, coef = list(mu, sig), trunc=1, log=TRUE))
pl.par <- poilogMLE(x1, startVals = c(mu = mean(log(x1)) + log(0.5), sig = sd(log(x1))))$par
nvl1.mle = mle2(nvl1, start=as.list(pl.par))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(mu, sig) -sum(dtrunc("poilog", x2, coef = list(mu, sig), trunc=1, log=TRUE))
pl.par <- poilogMLE(x2, startVals = c(mu = mean(log(x2)) + log(0.5), sig = sd(log(x2))))$par
nvl2.mle = mle2(nvl2, start=as.list(pl.par), method="CG")
summary(nvl2.mle)                            
logLik(nvl2.mle)

###teste com funcao fitpoilog
(tmodel1 <- fitpoilog(x1, trunc = 1))
summary(tmodel1)
(tmodel1 <- fitpoilog(x1, trunc = 1))
summary(tmodel1)
attributes(tmodel1)
(tmodel2 <- fitpoilog(x2, trunc = 1, method="CG", control = list(maxit = 2000)))

##################################################################################
### Testes da fitlnorm
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(meanlog, sdlog) -sum(dlnorm(x1, meanlog, sdlog, log = TRUE))
nvl1.mle <- mle2(nvl1, start = list(meanlog = mean(log(x1)), sdlog = sd(log(x1))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(meanlog, sdlog) -sum(dlnorm(x2, meanlog, sdlog, log=TRUE))
nvl2.mle <- mle2(nvl2, start = list(meanlog = mean(log(x2)), sdlog = sd(log(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitlnorm
(model1 <- fitlnorm(x1))
attributes(model1)
(model2 <- fitlnorm(x2))
attributes(model2)

##truncada em zero -> igual a sem truncagem
x1 <- samp1               
nvl1 = function(meanlog, sdlog) -sum(dtrunc("lnorm", x1, coef = list(meanlog, sdlog), trunc = 0, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(meanlog = mean(log(x2)), sdlog = sd(log(x2))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(meanlog, sdlog) -sum(dtrunc("lnorm", x2, coef = list(meanlog, sdlog), trunc = 0, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(meanlog = mean(log(x2)), sdlog = sd(log(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitlnorm
(model1 <- fitlnorm(x1, trunc = 0))
attributes(model1)
(model2 <- fitlnorm(x2, trunc = 0))
attributes(model2)

##truncada em um
x1 <- samp1[samp1 != 1]               
nvl1 = function(meanlog, sdlog) -sum(dtrunc("lnorm", x1, coef = list(meanlog, sdlog), trunc = 1, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(meanlog = mean(log(x2)), sdlog = sd(log(x2))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(meanlog, sdlog) -sum(dtrunc("lnorm", x2, coef = list(meanlog, sdlog), trunc = 1, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(meanlog = mean(log(x2)), sdlog = sd(log(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitlnorm
(model1 <- fitlnorm(x1, trunc = 1))
attributes(model1)
(model2 <- fitlnorm(x2, trunc = 1))
attributes(model2)

##################################################################################
### Testes da fitgamma
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(shape, rate) -sum(dgamma(x1, shape, rate, log = TRUE))
#chute inicial
ka <- (mean(x1)/sd(x1))^2
theta <- var(x1)/mean(x1)
kahat <- function(k, dados){
  eq <- length(dados)*(log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
}
shape <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x1)$root
scale <- mean(x1)/ka
nvl1.mle <- mle2(nvl1, start = list(shape = shape, rate = 1/scale))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(shape, rate) -sum(dgamma(x2, shape, rate, log = TRUE))
#chute inicial
ka <- (mean(x2)/sd(x2))^2
theta <- var(x2)/mean(x2)
ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x2)$root
theta <- mean(x2)/ka
nvl2.mle <- mle2(nvl2, start = list(shape = ka, rate = 1/theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgamma
(model1 <- fitgamma(x1))
attributes(model1)
fitdistr(x1, "gamma")
(model2 <- fitgamma(x2))
attributes(model2)
fitdistr(x2, "gamma")

##truncada em zero -> igual sem truncagem
x1 <- samp1               
nvl1 = function(shape, rate) -sum(dtrunc("gamma", x1, coef = list(shape, rate), trunc = 0, log=TRUE))
#chute inicial
ka <- (mean(x1)/sd(x1))^2
theta <- var(x1)/mean(x1)
ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x1)$root
theta <- mean(x1)/ka
nvl1.mle <- mle2(nvl1, start = list(shape = ka, rate = 1/theta))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(shape, rate) -sum(dtrunc("gamma", x2, coef = list(shape, rate), trunc = 0, log=TRUE))
#chute inicial
ka <- (mean(x2)/sd(x2))^2
theta <- var(x2)/mean(x2)
ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x2)$root
theta <- mean(x2)/ka
nvl2.mle = mle2(nvl2, start =list(shape = ka, rate = 1/theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgamma
(model1 <- fitgamma(x1, trunc = 0))
attributes(model1)
(model2 <- fitgamma(x2, trunc = 0))
attributes(model2)

##truncada em um
x1 <- samp1[samp1 != 1]               
nvl1 = function(shape, rate) -sum(dtrunc("gamma", x1, coef = list(shape, rate), trunc = 1, log=TRUE))
#chute inicial
ka <- (mean(x1)/sd(x1))^2
theta <- var(x1)/mean(x1)
ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x1)$root
theta <- mean(x1)/ka
nvl1.mle <- mle2(nvl1, start = list(shape = ka, rate = 1/theta), method="SANN")
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(shape, rate) -sum(dtrunc("gamma", x2, coef = list(shape, rate), trunc = 1, log=TRUE))
#chute inicial
ka <- (mean(x2)/sd(x2))^2
theta <- var(x2)/mean(x2)
ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x2)$root
theta <- mean(x2)/ka
nvl2.mle = mle2(nvl2, start =list(shape = ka, rate = 1/theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgamma
(model1 <- fitgamma(x1, trunc = 1, method = "CG", control = list(maxit = 1000)))
attributes(model1)
(model2 <- fitgamma(x2, trunc = 1))
attributes(model2)

##################################################################################
### Testes da fitgeom
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(prob) -sum(dgeom(x1, prob, log = TRUE))
nvl1.mle <- mle2(nvl1, start = list(prob = 1/(1+ mean(x1))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(prob) -sum(dgeom(x2, prob, log=TRUE))
nvl2.mle <- mle2(nvl2, start = list(prob = 1/(1 + mean(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgeom
(model1 <- fitgeom(x1))
attributes(model1)
(model2 <- fitgeom(x2))
attributes(model2)

##truncada em zero
x1 <- samp1               
nvl1 = function(prob) -sum(dtrunc("geom", x1, coef = prob, trunc = 0, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(prob = 1/(1 + mean(x2))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(prob) -sum(dtrunc("geom", x2, coef = prob, trunc = 0, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(prob = 1/(1 + mean(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgeom
(model1 <- fitgeom(x1, trunc = 0))
attributes(model1)
(model2 <- fitgeom(x2, trunc = 0))
attributes(model2)

##truncada em um
x1 <- samp1[samp1 != 1]               
nvl1 = function(prob) -sum(dtrunc("geom", x1, coef = prob, trunc = 1, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(prob = 1/(1 + mean(x2))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(prob) -sum(dtrunc("geom", x2, coef = prob, trunc = 1, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(prob = 1/(1 + mean(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitgeom
(model1 <- fitgeom(x1, trunc = 1, start = 0.016905))
attributes(model1)
(model2 <- fitgeom(x2, trunc = 1))
attributes(model2)

##################################################################################
### Testes da fitnbinom
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(prob) -sum(dnbinom(x1, size = length(x1), prob, log = TRUE))
nvl1.mle <- mle2(nvl1, start = list(prob = length(x1)/(length(x1) + mean(x1))), method="Brent", upper = 1, lower = 0)
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(prob) -sum(dnbinom(x2, size = length(x2), prob, log=TRUE))
nvl2.mle <- mle2(nvl2, start = list(prob = length(x2)/(length(x2) + mean(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitnbinom
(model1 <- fitnbinom(x1))
attributes(model1)
(model2 <- fitnbinom(x2))
attributes(model2)

##truncada em um
x1 <- samp1[samp1 != 1]               
nvl1 = function(prob) -sum(dtrunc("nbinom", x1, coef = list(size = length(x1), prob), trunc = 1, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(prob = length(x1)/(length(x1) + mean(x1))))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(prob) -sum(dtrunc("nbinom", x2, coef = list(size = length(x2), prob), trunc = 1, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(prob = length(x2)/(length(x2) + mean(x2))))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitnbinom
(model1 <- fitnbinom(x1, trunc = 1))
summary(model1)
attributes(model1)
(model2 <- fitnbinom(x2, trunc = 1))
attributes(model2)
summary(model2)

##################################################################################
### Testes da fitpower
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(s) -sum(dpower(x1, s, log = TRUE))
nvl1.mle <- mle2(nvl1, start = list(s = 2), method="Brent", upper = 20, lower = 1)
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(s) -sum(dpower(x2, s, log=TRUE))
nvl2.mle <- mle2(nvl2, start = list(s = 2), method="Brent", upper = 20, lower = 1)
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpower
(model1 <- fitpower(x1))
attributes(model1)
(model2 <- fitpower(x2))
attributes(model2)

##truncada em zero
x1 <- samp1[samp1 != 1]               
nvl1 = function(s) -sum(dtrunc("power", x1, coef = s, trunc = 1, log=TRUE))
nvl1.mle <- mle2(nvl1, start = list(s = 2), method="Brent", upper = 20, lower = 1)
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(s) -sum(dtrunc("power", x2, coef = s, trunc = 1, log=TRUE))
nvl2.mle = mle2(nvl2, start = list(s = 2), method="Brent", upper = 20, lower = 1)
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitpower
(model1 <- fitpower(x1, trunc = 1))
summary(model1)
attributes(model1)
(model2 <- fitpower(x2, trunc = 1))
attributes(model2)
summary(model2)

##################################################################################
### Testes da fitweibull
#################################################################################
## sem truncagem
x1 <- samp1               
nvl1 = function(shape, scale) -sum(dweibull(x1, shape, scale, log = TRUE))
#chute inicial
ka <- 1
theta <- mean(x1)
for(i in 1:100){
  theta <- (sum(x1^ka)/length(x1))^(1/ka)
  ka <- length(x1)/(sum(x1^ka * log(x1)) - sum(log(x1))/theta)
}
nvl1.mle <- mle2(nvl1, start = list(shape = ka, scale = theta))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(shape, scale) -sum(dweibull(x2, shape, scale, log = TRUE))
#chute inicial
ka <- 1
theta <- mean(x2)
for(i in 1:100){
  theta <- (sum(x2^ka)/length(x2))^(1/ka)
  ka <- length(x2)/(sum(x2^ka * log(x2)) - sum(log(x2))/theta)
}
nvl2.mle <- mle2(nvl2, start = list(shape = ka, scale = theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitweibull
(model1 <- fitweibull(x1))
attributes(model1)
fitdistr(x1, "weibull")
(model2 <- fitweibull(x2))
attributes(model2)
fitdistr(x2, "weibull")

##truncada em zero -> igual sem truncagem
x1 <- samp1               
nvl1 = function(shape, scale) -sum(dtrunc("weibull", x1, coef = list(shape, scale), trunc = 0, log=TRUE))
ka <- 1
theta <- mean(x1)
for(i in 1:100){
  theta <- (sum(x1^ka)/length(x1))^(1/ka)
  ka <- length(x1)/(sum(x1^ka * log(x1)) - sum(log(x1))/theta)
}
nvl1.mle <- mle2(nvl1, start = list(shape = ka, scale = theta))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2
nvl2 = function(shape, scale) -sum(dtrunc("weibull", x2, coef = list(shape, scale), trunc = 0, log=TRUE))
#chute inicial
ka <- 1
theta <- mean(x2)
for(i in 1:100){
  theta <- (sum(x2^ka)/length(x2))^(1/ka)
  ka <- length(x2)/(sum(x2^ka * log(x2)) - sum(log(x2))/theta)
}
nvl2.mle <- mle2(nvl2, start = list(shape = ka, scale = theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitweibull
(model1 <- fitweibull(x1, trunc = 0))
attributes(model1)
(model2 <- fitweibull(x2, trunc = 0))
attributes(model2)

##truncada em zero -> igual sem truncagem
x1 <- samp1[samp1 != 1]               
nvl1 = function(shape, scale) -sum(dtrunc("weibull", x1, coef = list(shape, scale), trunc = 1, log=TRUE))
ka <- 1
theta <- mean(x1)
for(i in 1:100){
  theta <- (sum(x1^ka)/length(x1))^(1/ka)
  ka <- length(x1)/(sum(x1^ka * log(x1)) - sum(log(x1))/theta)
}
nvl1.mle <- mle2(nvl1, start = list(shape = ka, scale = theta))
summary(nvl1.mle)                            
logLik(nvl1.mle)

x2 <- samp2[samp2 != 1]
nvl2 = function(shape, scale) -sum(dtrunc("weibull", x2, coef = list(shape, scale), trunc = 1, log=TRUE))
#chute inicial
ka <- 1
theta <- mean(x2)
for(i in 1:100){
  theta <- (sum(x2^ka)/length(x2))^(1/ka)
  ka <- length(x2)/(sum(x2^ka * log(x2)) - sum(log(x2))/theta)
}
nvl2.mle <- mle2(nvl2, start = list(shape = ka, scale = theta))
summary(nvl2.mle)
logLik(nvl2.mle)

###teste com funcao fitweibull
(model1 <- fitweibull(x1, trunc = 1))
attributes(model1)
(model2 <- fitweibull(x2, trunc = 1))
attributes(model2)