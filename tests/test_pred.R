library(sads)
library(bbmle)
library(poilog)
library(untb)

## Ajustes e graficos passo a passo
## Uma mostra Poisson de uma lognormal
samp1 <- rsad(100,frac=0.15,sad=lnorm,samp="Poisson",meanlog=3,sdlog=2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N=sum(samp1),S=length(samp1),nmax=sum(samp1))

## Ajustes dos modelos de sad

## Amostra 1
## Poilog
samp1.pln <- fitpoilog(samp1)
cf.pln <- as.numeric(samp1.pln@coef)
## Logserie
samp1.ls <- fitls(samp1)
cf.ls <- as.numeric(samp1.ls@coef)

## Amostra 2
## Poilog
samp2.pln <- fitpoilog(samp2)
cf2.pln <- as.numeric(samp2.pln@coef)
## ajuste à logserie
samp2.ls <- fitls(samp2)
cf2.ls <- as.numeric(samp2.ls@coef)


## Selecao de modelos
AICtab(samp1.pln, samp1.ls, nobs=length(samp1), base=T, weights=T)
AICtab(samp2.pln, samp2.ls, nobs=length(samp2), base=T, weights=T)## O AICctab nao funfa, por algum problema com o nobs: VERIFICAR


##GRAFICOS##
radpred2 <- function(object,...){
  psad <- paste("p", object@sad, sep="")
  dots <- c(as.list(object@coef), list(...))
  if(object@sad=="ls" & !"N" %in% names(dots)) dots$N <- sum(object@data$x)
  S <- length(object@data$x)
  y <- 1:max(object@data$x)
  X <- do.call(psad, c(list(q=y, lower.tail=F), dots))
  f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
  ab <- f1(ppoints(S))
  if(is.na(ab[1]) & !any(is.na(ab[-1]))){
    ab[1] <- sum(object@data$x) - sum(ab[-1])
  }
  new("rad", data.frame(rank=1:S, abund=ab))
}

## Amostra de uma lognormal##
## Esperados para as rank-plot com lista de coeficientes
samp1.pl.rad <- radpred(x=samp1, sad="poilog", coef=as.list(samp1.pln@coef))
samp1.ls.rad <- radpred(x=samp1, sad="ls", coef=as.list(samp1.ls@coef))
plot(rad(samp1))
points(samp1.pl.rad2)
points(samp1.ls.rad2, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

samp1.pl.rad2 <- radpred2(samp1.pln)
samp1.ls.rad2 <- radpred2(samp1.ls)
plot(rad(samp1))
points(samp1.pl.rad2)
points(samp1.ls.rad2, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))
#undebug(radpred2)

## Esperados para as oitavas com lista de coeficientes
octavpred2 <- function(object, oct, ...){
  S <- length(object@data$x)
  if(missing(oct)) oct <- 1:(ceiling(max(log2((object@data$x))))+1)
  n <- 2^(oct-1)
  psad <- paste("p", object@sad, sep = "")
  dots <- c(as.list(object@coef), list(...))
  if(object@sad == "ls" & !"N" %in% names(dots)) dots$N <- sum(object@data$x)
  Y <- do.call(psad, c(list(q = n), dots))
  Y <- c(Y[1], diff(Y))*S
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}

samp1.pl.octav <- octavpred(x=samp1, sad="poilog", coef=as.list(samp1.pln@coef))
samp1.ls.octav <- octavpred(x=samp1, sad="ls", coef=as.list(samp1.ls@coef))
plot(octav(samp1))
points(samp1.pl.octav)
points(samp1.ls.octav, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

samp1.pl.octav2 <- octavpred2(samp1.pln)
samp1.ls.octav2 <- octavpred2(samp1.ls)
plot(octav(samp1))
points(samp1.pl.octav2)
points(samp1.ls.octav2, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))
#debug(octavpred2)

## Amostra de uma logserie##
## Esperados para as rank-plot com lista de coeficientes
samp2.pl.rad <- radpred(x=samp2,sad="poilog",coef=as.list(samp2.pln@coef))
samp2.ls.rad <- radpred(x=samp2,sad="ls",coef=as.list(samp2.ls@coef))
plot(rad(samp2))
points(samp2.pl.rad)
points(samp2.ls.rad,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

samp2.pl.rad2 <- radpred2(samp2.pln)
samp2.ls.rad2 <- radpred2(samp2.ls)
plot(rad(samp2))
points(samp2.pl.rad2)
points(samp2.ls.rad2, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

## Esperados para as oitavas com lista de coeficientes
samp2.pl.octav <- octavpred(x=samp2,sad="poilog",coef=as.list(samp2.pln@coef))
samp2.ls.octav <- octavpred(x=samp2,sad="ls",coef=as.list(samp2.ls@coef))
plot(octav(samp2),ylim=c(0,max(c(octav(samp2)$Freq, samp2.pl.octav$Freq,samp2.ls.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))
                            
samp2.pl.octav2 <- octavpred2(samp2.pln)
samp2.ls.octav2 <- octavpred2(samp2.ls)
plot(octav(samp2),ylim=c(0,max(c(octav(samp2)$Freq, samp2.pl.octav$Freq,samp2.ls.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

##Amostra biomass
aa1 <- runif(50, 2^-6, 2^-5)
aa2 <- runif(25, 2^-5, 2^-4)
aa3 <- runif(13, 2^-4, 2^-3)
aa4 <- runif(7, 2^-3, 2^-2)
aa5 <- runif(4, 2^-2, 2^-1)
aa6 <- runif(1, 2^-1, 2^0)
aa <- c(aa1, aa2, aa3, aa4, aa5, aa6)

## Amostra 1
## Poilog
aa.ln <- fitlnorm(aa)
cf.ln <- as.numeric(aa.ln@coef)
## Logserie
aa.bu <- fitweibull(aa)
cf.bu <- as.numeric(aa.bu@coef)

AICtab(aa.ln, aa.bu, nobs=length(aa), base=T, weights=T)

aa.ln.rad2 <- radpred2(aa.ln)
aa.bu.rad2 <- radpred2(aa.bu)
plot(rad(aa))
plot(aa.ln.rad2)
points(aa.bu.rad2, col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))