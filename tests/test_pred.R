library(sads)
#library(bbmle)
#library(poilog)
library(untb)

## Ajustes e graficos passo a passo
## Uma mostra Poisson de uma lognormal
set.seed(1913)
samp1 <- rsad(200, frac=0.15, sad=lnorm, samp="Poisson", meanlog=3, sdlog=2)
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N=sum(samp1), S=length(samp1), nmax=sum(samp1))

## Ajustes dos modelos de sad

## Amostra 1
## Poilog
samp1.pln <- fitpoilog(samp1,trunc=0)
## Logserie
samp1.ls <- fitls(samp1)
## Gamma
samp1.gm <- fitgamma(samp1, trunc = 0)
## Power
samp1.pw <- fitpower(samp1)

## Amostra 2
## Poilog
samp2.pln <- fitpoilog(samp2, trunc=0)
## ajuste à logserie
samp2.ls <- fitls(samp2)
## Gamma
samp2.gm <- fitgamma(samp2, trunc=0.5)
## Power
samp2.pw <- fitpower(samp2)


## Selecao de modelos
AICtab(samp1.pln, samp1.ls, samp1.gm, samp1.pw, nobs=length(samp1), base=T, weights=T)
AICtab(samp2.pln, samp2.ls, samp2.gm, samp2.pw, nobs=length(samp2), base=T, weights=T)

##GRAFICOS##
radpred2 <- function(object,...){
  dots <- list(...)
  S <- length(object@data$x)
  y <- 1:max(object@data$x)
  if(!is.na(object@trunc)){
    if(object@sad=="ls") 
      X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(alpha = object@coef, N = sum(object@data$x)), lower.tail=F, trunc = object@trunc), dots))
    else
      X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
  }else {
    psad <- get(paste("p", object@sad, sep=""), mode = "function")
    if(object@sad=="ls")
      X <- do.call(psad, c(list(q = y, lower.tail = F, alpha = object@coef, N = sum(object@data$x)), dots))
    else
      X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(object@coef), dots))
  }
  f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
  ab <- f1(ppoints(S))
  if(is.na(ab[1]) & !any(is.na(ab[-1]))){
    ab[1] <- sum(object@data$x) - sum(ab[-1])
  }
  new("rad", data.frame(rank=1:S, abund=ab))
}

## Nova versao radpred2 com funcao de quantil: muuuito lenta
radpred2 <- function(object,...){
  dots <- list(...)
  S <- length(object@data$x)
  y <- ppoints(S)
  if(!is.na(object@trunc)){
    if(object@sad=="ls") 
      X <- do.call(qtrunc, c(list(object@sad, p = y, coef = list(alpha = object@coef, N = sum(object@data$x)), lower.tail=F, trunc = object@trunc), dots))
    else
      X <- do.call(qtrunc, c(list(object@sad, p = y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
  }else {
    qsad <- get(paste("q", object@sad, sep=""), mode = "function")
    if(object@sad=="ls")
      X <- do.call(qsad, c(list(p = y, lower.tail = F, alpha = object@coef, N = sum(object@data$x)), dots))
    else
      X <- do.call(qsad, c(list(p = y, lower.tail = F), as.list(object@coef), dots))
  }
  
##  if(is.na(ab[1]) & !any(is.na(ab[-1]))){
##    ab[1] <- sum(object@data$x) - sum(ab[-1])
##  }
  #X[1] <- 
  new("rad", data.frame(rank=1:S, abund=X))
}

## Amostra de uma lognormal##
## Esperados para as rank-plot com lista de coeficientes
samp1.pl.rad <- radpred(x=samp1, sad="poilog", coef=as.list(coef(samp1.pln)))
samp1.ls.rad <- radpred(x=samp1, sad="ls", coef=as.list(coef(samp1.ls)))
samp1.gm.rad <- radpred(x=samp1, sad="gamma", coef=as.list(coef(samp1.gm)))
samp1.pw.rad <- radpred(x=samp1, sad="power", coef=as.list(samp1.pw@coef))
plot(rad(samp1), ylim=c(1, 500), main="Rad1 - samp1")
points(samp1.pl.rad)
points(samp1.ls.rad, col="red")
points(samp1.gm.rad, col="green")
points(samp1.pw.rad, col="purple")
legend("topright",  c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

samp1.pl.rad2 <- radpred2(samp1.pln)
samp1.ls.rad2 <- radpred2(samp1.ls)
samp1.gm.rad2 <- radpred2(samp1.gm)
samp1.pw.rad2 <- radpred2(samp1.pw)
plot(rad(samp1), ylim=c(1, 500), main="Rad2 - samp1")
points(samp1.pl.rad2)
points(samp1.ls.rad2, col="red")
points(samp1.gm.rad2, col="green")
points(samp1.pw.rad2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))
#undebug(radpred2)

## Sample 2
samp2.pl.rad2 <- radpred2(samp2.pln)
samp2.ls.rad2 <- radpred2(samp2.ls)
samp2.gm.rad2 <- radpred2(samp2.gm)
samp2.pw.rad2 <- radpred2(samp2.pw)
plot(rad(samp2), ylim=c(1, 500), main="Rad2 - samp2")
points(samp2.pl.rad2)
points(samp2.ls.rad2, col="red")
points(samp2.gm.rad2, col="green")
points(samp2.pw.rad2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))
## Estranho: a log-serie tem AIC similar a gamma, mas o fit parece mto pior
AIC(samp2.ls)
AIC(samp2.gm)
## Conferindo quantis
(S <- length(samp2))
Sp <- ppoints(S)
## truncagem em um
x1 <- samp1[samp1 != 1]
x2 <- samp2[samp2 != 1]
## Amostra 1
## Poilog
x1.pln <- fitpoilog(x1, trunc = 1)
## Logserie
x1.ls <- fitls(x1, trunc = 1)
## Gamma
x1.gm <- fitgamma(x1, trunc = 1)
## Power
x1.pw <- fitpower(x1, trunc = 1)

## Amostra 2
## Poilog
x2.pln <- fitpoilog(x2, trunc = 1)
## ajuste à logserie
x2.ls <- fitls(x2, trunc = 1)
## Gamma
x2.gm <- fitgamma(x2, trunc = 1)
## Power
x2.pw <- fitpower(x2, trunc = 1)

AICtab(x1.pln, x1.ls, x1.gm, x1.pw, nobs=length(x1), base=T, weights=T)
AICtab(x2.pln, x2.ls, x2.gm, x2.pw, nobs=length(x2), base=T, weights=T)

x1.pl.rad <- radpred(x=x1, sad="poilog", coef=as.list(x1.pln@coef))
x1.ls.rad <- radpred(x=x1, sad="ls", coef=as.list(x1.ls@coef))
x1.gm.rad <- radpred(x=x1, sad="gamma", coef=as.list(x1.gm@coef))
x1.pw.rad <- radpred(x=x1, sad="power", coef=as.list(x1.pw@coef))
plot(rad(samp1), ylim=c(1, 500), main="Rad1 - x1")
points(x1.pl.rad)
points(x1.ls.rad, col="red")
points(x1.gm.rad, col="green")
points(x1.pw.rad, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

x1.pl.rad2 <- radpred2(x1.pln)
x1.ls.rad2 <- radpred2(x1.ls)
x1.gm.rad2 <- radpred2(x1.gm)
x1.pw.rad2 <- radpred2(x1.pw)
plot(rad(samp1), ylim=c(1, 500), main="Rad2 - x1")
points(x1.pl.rad2)
points(x1.ls.rad2, col="red")
points(x1.gm.rad2, col="green")
points(x1.pw.rad2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

## Esperados para as oitavas com lista de coeficientes
octavpred2 <- function(object, oct, ...){
  dots <- c(list(...))
  S <- length(object@data$x)
  if(missing(oct)) oct <- 1:(ceiling(max(log2((object@data$x))))+1)
  n <- 2^(oct-1)
  if(!is.na(object@trunc)){
    if(object@sad == "ls")
      Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = list(alpha = object@coef, N = sum(object@data$x)), trunc = object@trunc), dots))
    else
      Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = as.list(object@coef), trunc = object@trunc), dots))
  }else {
    psad <- get(paste("p", object@sad, sep=""), mode = "function")
    if(object@sad == "ls")
      Y <- do.call(psad, c(list(q = n, alpha = object@coef, N = sum(object@data$x)), dots))
    else
      Y <- do.call(psad, c(list(q = n), as.list(object@coef), dots))
  }
  Y <- c(Y[1], diff(Y))*S
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}

samp1.pl.octav <- octavpred(x=samp1, sad="poilog", coef=as.list(samp1.pln@coef))
samp1.ls.octav <- octavpred(x=samp1, sad="ls", coef=as.list(samp1.ls@coef))
samp1.gm.octav <- octavpred(x=samp1, sad="gamma", coef=as.list(samp1.gm@coef))
samp1.pw.octav <- octavpred(x=samp1, sad="power", coef=as.list(samp1.pw@coef))
plot(octav(samp1), ylim=c(0, 25))
points(samp1.pl.octav)
points(samp1.ls.octav, col="red")
points(samp1.gm.octav, col="green")
points(samp1.pw.octav, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

samp1.pl.octav2 <- octavpred2(samp1.pln)
samp1.ls.octav2 <- octavpred2(samp1.ls)
samp1.gm.octav2 <- octavpred2(samp1.gm)
samp1.pw.octav2 <- octavpred2(samp1.pw)
plot(octav(samp1), ylim=c(0, 25))
points(samp1.pl.octav2)
points(samp1.ls.octav2, col="red")
points(samp1.gm.octav2, col="green")
points(samp1.pw.octav2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))
#debug(octavpred2)

## truncagem em um
x1.pl.octav <- octavpred(x=x1, sad="poilog", coef=as.list(x1.pln@coef))
x1.ls.octav <- octavpred(x=x1, sad="ls", coef=as.list(x1.ls@coef))
x1.gm.octav <- octavpred(x=x1, sad="gamma", coef=as.list(x1.gm@coef))
x1.pw.octav <- octavpred(x=x1, sad="power", coef=as.list(x1.pw@coef))
plot(octav(x1), ylim=c(0, 35))
points(x1.pl.octav)
points(x1.ls.octav, col="red")
points(x1.gm.octav, col="green")
points(x1.pw.octav, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

x1.pl.octav2 <- octavpred2(x1.pln)
x1.ls.octav2 <- octavpred2(x1.ls)
x1.gm.octav2 <- octavpred2(x1.gm)
x1.pw.octav2 <- octavpred2(x1.pw)
plot(octav(x1), ylim=c(0, 35))
points(x1.pl.octav2)
points(x1.ls.octav2, col="red")
points(x1.gm.octav2, col="green")
points(x1.pw.octav2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

## Amostra de uma logserie##
## Esperados para as rank-plot com lista de coeficientes
samp2.pl.rad <- radpred(x=samp2,sad="poilog",coef=as.list(samp2.pln@coef))
samp2.ls.rad <- radpred(x=samp2,sad="ls",coef=as.list(samp2.ls@coef))
samp2.gm.rad <- radpred(x=samp2,sad="gamma",coef=as.list(samp2.gm@coef))
samp2.pw.rad <- radpred(x=samp2,sad="power",coef=as.list(samp2.pw@coef))
plot(rad(samp2), ylim=c(1, 500), main= "Rad1 - samp2")
points(samp2.pl.rad)
points(samp2.ls.rad,col="red")
points(samp2.gm.rad,col="green")
points(samp2.pw.rad,col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

samp2.pl.rad2 <- radpred2(samp2.pln)
samp2.ls.rad2 <- radpred2(samp2.ls)
samp2.gm.rad2 <- radpred2(samp2.gm)
samp2.pw.rad2 <- radpred2(samp2.pw)
plot(rad(samp2), ylim=c(1, 500), main= "Rad2 - samp2")
points(samp2.pl.rad2)
points(samp2.ls.rad2, col="red")
points(samp2.gm.rad2, col="green")
points(samp2.pw.rad2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

##Truncados em um
x2.pl.rad <- radpred(x=x2, sad="poilog", coef=as.list(x2.pln@coef))
x2.ls.rad <- radpred(x=x2, sad="ls", coef=as.list(x2.ls@coef))
x2.gm.rad <- radpred(x=x2, sad="gamma", coef=as.list(x2.gm@coef))
x2.pw.rad <- radpred(x=x2, sad="power", coef=as.list(x2.pw@coef))
plot(rad(samp2), ylim=c(1, 700), main= "Rad1 - x2")
points(x2.pl.rad)
points(x2.ls.rad, col="red")
points(x2.gm.rad, col="green")
points(x2.pw.rad, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

x2.pl.rad2 <- radpred2(x2.pln)
x2.ls.rad2 <- radpred2(x2.ls)
x2.gm.rad2 <- radpred2(x2.gm)
x2.pw.rad2 <- radpred2(x2.pw)
plot(rad(samp2), ylim=c(1, 700), main= "Rad2 - x2")
points(x2.pl.rad2)
points(x2.ls.rad2, col="red")
points(x2.gm.rad2, col="green")
points(x2.pw.rad2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

## Esperados para as oitavas com lista de coeficientes
samp2.pl.octav <- octavpred(x=samp2, sad="poilog", coef=as.list(samp2.pln@coef))
samp2.ls.octav <- octavpred(x=samp2, sad="ls", coef=as.list(samp2.ls@coef))
samp2.gm.octav <- octavpred(x=samp2, sad="gamma", coef=as.list(samp2.gm@coef))
samp2.pw.octav <- octavpred(x=samp2, sad="power", coef=as.list(samp2.pw@coef))

plot(octav(samp2), main="Oct1 - samp2", ylim=c(0, 5+max(c(octav(samp2)$Freq, samp2.pl.octav$Freq, samp2.ls.octav$Freq, samp2.gm.octav$Freq, samp2.pw.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav, col = "red")
points(samp2.gm.octav, col = "green")
points(samp2.pw.octav, col = "purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))
                            
samp2.pl.octav2 <- octavpred2(samp2.pln)
samp2.ls.octav2 <- octavpred2(samp2.ls)
samp2.gm.octav2 <- octavpred2(samp2.gm)
samp2.pw.octav2 <- octavpred2(samp2.pw)
plot(octav(samp2), main="Oct2 - samp2", ylim=c(0, 5+max(c(octav(samp2)$Freq, samp2.pl.octav$Freq, samp2.ls.octav$Freq, samp2.gm.octav$Freq, samp2.pw.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav,col="red")
points(samp2.gm.octav,col="green")
points(samp2.pw.octav,col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

## truncagem em um
x2.pl.octav <- octavpred(x=x2, sad="poilog", coef=as.list(x2.pln@coef))
x2.ls.octav <- octavpred(x=x2, sad="ls", coef=as.list(x2.ls@coef))
x2.gm.octav <- octavpred(x=x2, sad="gamma", coef=as.list(x2.gm@coef))
x2.pw.octav <- octavpred(x=x2, sad="power", coef=as.list(x2.pw@coef))
plot(octav(x2), main="Oct1 - x2", ylim =c(0, 25))
points(x2.pl.octav)
points(x2.ls.octav, col="red")
points(x2.gm.octav, col="green")
points(x2.pw.octav, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

x2.pl.octav2 <- octavpred2(x2.pln)
x2.ls.octav2 <- octavpred2(x2.ls)
x2.gm.octav2 <- octavpred2(x2.gm)
x2.pw.octav2 <- octavpred2(x2.pw)
plot(octav(x2), main = "Oct2 - x2", ylim=c(0, 25))
points(x2.pl.octav2)
points(x2.ls.octav2, col="red")
points(x2.gm.octav2, col="green")
points(x2.pw.octav2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

##Amostra biomass
radpredc2 <- function(object, ...){
  dots <- list(...)
  S <- length(object@data$x)
  Y <- ppoints(S)
  if(!is.na(object@trunc)){
    ab <- do.call(qtrunc, c(list(object@sad, p = Y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
  }else{
    qsad <- get(paste("q", object@sad, sep=""), mode = "function")
    ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(object@coef), dots))
  }
  new("rad",data.frame(rank=1:S, abund=ab))
}

en.b <- read.table("en_b.txt")
en.b <- en.b$V2
bes.ln <- fitlnorm(en.b)
bes.gm <- fitgamma(en.b, star = c(0.4471, 0.6119))
bes.we <- fitweibull(en.b)


bes.ln.rad2 <- radpredc2(bes.ln)
bes.gm.rad2 <- radpredc2(bes.gm)
bes.we.rad2 <- radpredc2(bes.we)
plot(rad(en.b), main="Rad cont",ylim=c(min(en.b, bes.ln.rad2$abund, bes.gm.rad2$abund, bes.we.rad2$abund), max(en.b, bes.ln.rad2$abund, bes.gm.rad2$abund, bes.we.rad2$abund)))
points(bes.ln.rad2)
points(bes.gm.rad2, col="red")
points(bes.we.rad2, col="green")
legend("topright", c("Logseries", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green", "purple"))

en.bt <- en.b[en.b > 0.05]
bes.lnt <- fitlnorm(en.bt, trunc = 0.05)
bes.gmt <- fitgamma(en.bt, star = c(0.56, 0.61), trunc=0.05)
bes.wet <- fitweibull(en.bt, star=c(0.79, 0.87), trunc =0.05)

bes.ln.rad2t <- radpredc2(bes.lnt)
bes.gm.rad2t <- radpredc2(bes.gmt)
bes.we.rad2t <- radpredc2(bes.wet)
plot(rad(en.b), main="Rad2 cont",ylim=c(min(en.b, bes.ln.rad2t$abund, bes.gm.rad2t$abund, bes.we.rad2t$abund), max(en.b, bes.ln.rad2t$abund, bes.gm.rad2t$abund, bes.we.rad2t$abund)))
points(bes.ln.rad2t)
points(bes.gm.rad2t, col="red")
points(bes.we.rad2t, col="green")
legend("topright", c("Logseries", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green", "purple"))

radpredt <- function(object, ...){
  dots <- list(...)
  S <- length(object@data$x)
  if (object@sad == "ls" || object@sad == "geom" || object@sad == "nbinom"|| object@sad == "zipf"|| object@sad == "power"|| object@sad == "poilog"){
    y <- 1:max(object@data$x)
    if(!is.na(object@trunc)){
      if(object@sad=="ls") 
        X <- do.call(ptrunc, c(list(object@sad, q = y, coef = list(alpha = object@coef, N = sum(object@data$x)), lower.tail=F, trunc = object@trunc), dots))
      else
        X <- do.call(ptrunc, c(list(object@sad, q = y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
    }else {
      psad <- get(paste("p", object@sad, sep=""), mode = "function")
      if(object@sad=="ls")
        X <- do.call(psad, c(list(q = y, lower.tail = F, alpha = object@coef, N = sum(object@data$x)), dots))
      else
        X <- do.call(psad, c(list(q = y, lower.tail = F), as.list(object@coef), dots))
    }
    f1 <- approxfun(x=c(1, X), y=c(0, y), method="constant")
    ab <- f1(ppoints(S))
    if(is.na(ab[1]) & !any(is.na(ab[-1]))){
      ab[1] <- sum(object@data$x) - sum(ab[-1])
    }
  }else if(object@sad == "gamma" || object@sad == "lnorm" || object@sad == "weibull"){
    Y <- ppoints(S)
    if(!is.na(object@trunc)){
      ab <- do.call(qtrunc, c(list(object@sad, p = Y, coef = as.list(object@coef), lower.tail=F, trunc = object@trunc), dots))
    }else{
      qsad <- get(paste("q", object@sad, sep=""), mode = "function")
      ab <- do.call(qsad, c(list(p=Y, lower.tail=F), as.list(object@coef), dots))
    }
  } else
     stop("unsupported distribution")
  new("rad", data.frame(rank=1:S, abund=ab))
}
#CONTINUA
##sem truncagem
bes.ln.radt <- radpredt(bes.ln)
bes.gm.radt <- radpredt(bes.gm)
bes.we.radt <- radpredt(bes.we)
plot(rad(en.b), main="Rad t", ylim=c(min(en.b, bes.ln.radt$abund, bes.gm.radt$abund, bes.we.radt$abund), max(en.b, bes.ln.radt$abund, bes.gm.radt$abund, bes.we.radt$abund)))
points(bes.ln.radt)
points(bes.gm.radt, col="red")
points(bes.we.radt, col="green")
legend("topright", c("Logseries", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green", "purple"))

##com truncagem
bes.ln.radt <- radpredt(bes.lnt)
bes.gm.radt <- radpredt(bes.gmt)
bes.we.radt <- radpredt(bes.wet)
plot(rad(en.b), main="Rad t trunc",ylim=c(min(en.b, bes.ln.radt$abund, bes.gm.radt$abund, bes.we.radt$abund), max(en.b, bes.ln.radt$abund, bes.gm.radt$abund, bes.we.radt$abund)))
points(bes.ln.radt)
points(bes.gm.radt, col="red")
points(bes.we.radt, col="green")
legend("topright", c("Logseries", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green", "purple"))

#DISCRETA
##sem truncagem
samp1.pl.radt <- radpredt(samp1.pln)
samp1.ls.radt <- radpredt(samp1.ls)
samp1.pw.radt <- radpredt(samp1.pw)
plot(rad(samp1), ylim=c(1, 500), main="Radt - samp1")
points(samp1.pl.radt)
points(samp1.ls.radt, col="red")
points(samp1.pw.radt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Power"), lty=1, col=c("blue","red", "purple"))

##com truncagem
x1.pl.radt <- radpredt(x1.pln)
x1.ls.radt <- radpredt(x1.ls)
x1.pw.radt <- radpredt(x1.pw)
plot(rad(samp1), ylim=c(1, 500), main="Radt - x1")
points(x1.pl.radt)
points(x1.ls.radt, col="red")
points(x1.pw.radt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Power"), lty=1, col=c("blue","red", "purple"))

#dada pelo usuario
us.pl.radt <- radpredt(x = samp1, coef=c(0.8782991, 2.4773384), sad = "poilog")
us.ls.radt <- radpredt(x = samp1, coef=c(sum(samp1), 14.86564), sad = "ls")
us.pw.radt <- radpredt(x = samp1, coef=c(1.379446), sad = "power")
plot(rad(samp1), ylim=c(1, 500), main="Radt - samp1")
points(us.pl.radt)
points(us.ls.radt, col="red")
points(us.pw.radt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Power"), lty=1, col=c("blue","red", "purple"))


octavpred2 <- function(object, oct, ...){
  dots <- list(...)
  S <- length(object@data$x)
  if(missing(oct)){
    oct <- 1:(ceiling(max(log2(object@data$x)))+1)
    if(any(object@data$x < 1)){
      octlower <- ceiling(min(log2((object@data$x)))+1):0
      oct <- c(octlower, oct)
    }
  }
  n <- 2^(oct-1)
  if(!is.na(object@trunc)){
    if(object@sad == "ls")
      Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = list(alpha = object@coef, N = sum(object@data$x)), trunc = object@trunc), dots))
    else
      Y <- do.call(ptrunc, c(list(object@sad, q = n, coef = as.list(object@coef), trunc = object@trunc), dots))
  }else {
    psad <- get(paste("p", object@sad, sep=""), mode = "function")
    if(object@sad == "ls")
      Y <- do.call(psad, c(list(q = n, alpha = object@coef, N = sum(object@data$x)), dots))
    else
      Y <- do.call(psad, c(list(q = n), as.list(object@coef), dots))
  }
  Y <- c(Y[1], diff(Y))*S
  new("octav", data.frame(octave = oct, upper = factor(n), Freq = Y))
}

#DISCRETA
#sem truncagem
samp1.pl.octav2 <- octavpred2(samp1.pln)
samp1.ls.octav2 <- octavpred2(samp1.ls)
samp1.gm.octav2 <- octavpred2(samp1.gm)
samp1.pw.octav2 <- octavpred2(samp1.pw)
plot(octav(samp1), ylim=c(0, 25))
points(samp1.pl.octav2)
points(samp1.ls.octav2, col="red")
points(samp1.gm.octav2, col="green")
points(samp1.pw.octav2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

#com truncagem
x1.pl.octav2 <- octavpred2(x1.pln)
x1.ls.octav2 <- octavpred2(x1.ls)
x1.gm.octav2 <- octavpred2(x1.gm)
x1.pw.octav2 <- octavpred2(x1.pw)
plot(octav(x1), ylim=c(0, 25))
points(x1.pl.octav2)
points(x1.ls.octav2, col="red")
points(x1.gm.octav2, col="green")
points(x1.pw.octav2, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

#CONTINUA
#sem truncagem
bes.ln.octav2 <- octavpred2(bes.ln)
bes.gm.octav2 <- octavpred2(bes.gm)
bes.we.octav2 <- octavpred2(bes.we)
bes.oct <- octav(en.b)
X <- c((min(as.integer(as.character(bes.oct$octave)))-1), as.integer(as.character(bes.oct$octave)))
X <- X[-length(X)]+diff(X)/2
plot(X, bes.oct$Freq, type="h", ylim=c(0, max(bes.octt$Freq)) , xlab="Abundance class", ylab="Frequency")
points(bes.ln.octav2)
points(bes.gm.octav2, col="red")
points(bes.we.octav2, col = "green")
legend("topleft", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

bes.hist<-rep(as.integer(as.character(bes.oct$octave)), as.integer(as.character(bes.oct$Freq)))
hist(bes.hist, col="gray", main="", ylab= "Frequency", xlab="Abundance class", breaks=c(min(as.integer(as.character(bes.oct$octave)))-1, as.integer(as.character(bes.oct$octave))))

points(bes.ln.octav2)
points(bes.gm.octav2, col="red")
points(bes.we.octav2, col = "green")
legend("topleft", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

plot(octav(en.b))
plot(octav(en.b), ylim=c(0, 5+max(c(octav(en.b)$Freq, bes.ln.octav2$Freq, bes.gm.octav2$Freq, bes.we.octav2$Freq))))

#com truncagem
bes.ln.octav2t <- octavpred2(bes.lnt)
bes.gm.octav2t <- octavpred2(bes.gmt)
bes.we.octav2t <- octavpred2(bes.wet)
bes.octt <- octav(en.bt)
Xt <- c((min(as.integer(as.character(bes.octt$octave)))-1), as.integer(as.character(bes.octt$octave)))
Xt <- Xt[-length(Xt)]+diff(Xt)/2
plot(Xt, bes.octt$Freq, type="h", ylim=c(0, max(bes.octt$Freq)), xlab="Abundance class", ylab="Frequency")
points(bes.ln.octav2t)
points(bes.gm.octav2t, col="red")
points(bes.we.octav2t, col = "green")
legend("topright", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

bes.histt<-rep(as.integer(as.character(bes.octt$octave)), as.integer(as.character(bes.octt$Freq)))
hist(bes.histt, col="gray", breaks=c((min(as.integer(as.character(bes.octt$octave)))-1), as.integer(as.character(bes.octt$octave))))
points(bes.ln.octav2t)
points(bes.gm.octav2t, col="red")
points(bes.we.octav2t, col = "green")
legend("topright", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

plot(octav(en.bt))
plot(octav(en.bt), ylim=c(0, max(c(octav(en.bt)$Freq, bes.ln.octav2t$Freq, bes.gm.octav2t$Freq, bes.we.octav2t$Freq))))

#OCTAVPREDT
#DISCRETA
#sem truncagem
samp1.pl.octavt <- octavpredt(samp1.pln)
samp1.ls.octavt <- octavpredt(samp1.ls)
samp1.gm.octavt <- octavpredt(samp1.gm)
samp1.pw.octavt <- octavpredt(samp1.pw)
plot(octav(samp1), ylim=c(0, 25))
points(samp1.pl.octavt)
points(samp1.ls.octavt, col="red")
points(samp1.gm.octavt, col="green")
points(samp1.pw.octavt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

us.pl.octavt <- octavpredt(x = samp1, coef =c(0.8782991, 2.4773384), sad ="poilog", trunc = 0)
us.ls.octavt <- octavpredt(x = samp1, coef =c(sum(samp1), 14.86564), sad = "ls")
us.gm.octavt <- octavpredt(x = samp1, coef =c(0.47375785, 0.01467221), sad = "gamma", trunc = 0)
us.pw.octavt <- octavpredt(x = samp1, coef = 1.379446, sad = "power")
plot(octav(samp1), ylim=c(0, 25))
points(us.pl.octavt)
points(us.ls.octavt, col="red")
points(us.gm.octavt, col="green")
points(us.pw.octavt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

#com truncagem
x1.pl.octavt <- octavpredt(x1.pln)
x1.ls.octavt <- octavpredt(x1.ls)
x1.gm.octavt <- octavpredt(x1.gm)
x1.pw.octavt <- octavpredt(x1.pw)
plot(octav(x1), ylim=c(0, 25))
points(x1.pl.octavt)
points(x1.ls.octavt, col="red")
points(x1.gm.octavt, col="green")
points(x1.pw.octavt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))

x1.us.pl.octavt <- octavpredt(x = x1, coef =c(1.635917, 2.135059), sad ="poilog", trunc = 1)
x1.us.ls.octavt <- octavpredt(x = x1, coef =c(sum(x1), 14.39427), sad = "ls", trunc = 1)
x1.us.gm.octavt <- octavpredt(x = x1, coef =c(0.170643892, 0.007949086), sad = "gamma", trunc = 1)
x1.us.pw.octavt <- octavpredt(x = x1, coef = 1.444122, sad = "power", trunc =1)
plot(octav(x1), ylim=c(0, 25))
points(x1.us.pl.octavt)
points(x1.us.ls.octavt, col="red")
points(x1.us.gm.octavt, col="green")
points(x1.us.pw.octavt, col="purple")
legend("topright", c("Poisson-lognormal", "Logseries", "Gamma", "Power"), lty=1, col=c("blue","red", "green", "purple"))


#CONTINUA
#sem truncagem
bes.ln.octav2 <- octavpred2(bes.ln)
bes.gm.octav2 <- octavpred2(bes.gm)
bes.we.octav2 <- octavpred2(bes.we)
bes.oct <- octav(en.b)
X <- c((min(as.integer(as.character(bes.oct$octave)))-1), as.integer(as.character(bes.oct$octave)))
X <- X[-length(X)]+diff(X)/2
plot(X, bes.oct$Freq, type="h", ylim=c(0, max(bes.octt$Freq)) , xlab="Abundance class", ylab="Frequency")
points(bes.ln.octav2)
points(bes.gm.octav2, col="red")
points(bes.we.octav2, col = "green")
legend("topleft", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

bes.hist<-rep(as.integer(as.character(bes.oct$octave)), as.integer(as.character(bes.oct$Freq)))
hist(bes.hist, col="gray", main="", ylab= "Frequency", xlab="Abundance class", breaks=c(min(as.integer(as.character(bes.oct$octave)))-1, as.integer(as.character(bes.oct$octave))))

points(bes.ln.octav2)
points(bes.gm.octav2, col="red")
points(bes.we.octav2, col = "green")
legend("topleft", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

plot(octav(en.b))
plot(octav(en.b), ylim=c(0, 5+max(c(octav(en.b)$Freq, bes.ln.octav2$Freq, bes.gm.octav2$Freq, bes.we.octav2$Freq))))

bes.ln.octavt <- octavpredt(bes.ln, coef=c(-1.759404, 2.066396), sad = "lnorm")
bes.gm.octavt <- octavpredt(bes.gm, coef=c(0.4471375,0.6119233), sad= "gamma")
bes.we.octavt <- octavpredt(bes.we, coef=c(0.5711702,0.4705112), sad="weibull")
plot(octav(en.b))
points(bes.ln.octavt)
points(bes.gm.octavt, col="red")
points(bes.we.octavt, col = "green")
legend("topleft", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))


#com truncagem
bes.ln.octav2t <- octavpred2(bes.lnt)
bes.gm.octav2t <- octavpred2(bes.gmt)
bes.we.octav2t <- octavpred2(bes.wet)
bes.octt <- octav(en.bt)
Xt <- c((min(as.integer(as.character(bes.octt$octave)))-1), as.integer(as.character(bes.octt$octave)))
Xt <- Xt[-length(Xt)]+diff(Xt)/2
plot(Xt, bes.octt$Freq, type="h", ylim=c(0, max(bes.octt$Freq)), xlab="Abundance class", ylab="Frequency")
points(bes.ln.octav2t)
points(bes.gm.octav2t, col="red")
points(bes.we.octav2t, col = "green")
legend("topright", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

bes.histt<-rep(as.integer(as.character(bes.octt$octave)), as.integer(as.character(bes.octt$Freq)))
hist(bes.histt, col="gray", breaks=c((min(as.integer(as.character(bes.octt$octave)))-1), as.integer(as.character(bes.octt$octave))))
points(bes.ln.octav2t)
points(bes.gm.octav2t, col="red")
points(bes.we.octav2t, col = "green")
legend("topright", c("Lognormal", "Gamma", "Weibull"), lty=1, col=c("blue","red", "green"))

plot(octav(en.bt))
plot(octav(en.bt), ylim=c(0, max(c(octav(en.bt)$Freq, bes.ln.octav2t$Freq, bes.gm.octav2t$Freq, bes.we.octav2t$Freq))))
