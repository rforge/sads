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
cf.pln <- as.numeric(coef(samp1.pln))
## Logserie
samp1.ls <- fitlogser(samp1)
cf.ls <- as.numeric(coef(samp1.ls))

## Amostra 2
## Poilog
samp2.pln <- fitpoilog(samp2)
cf2.pln <- as.numeric(coef(samp2.pln))
## ajuste à logserie
samp2.ls <- fitlogser(samp2)
cf2.ls <- as.numeric(coef(samp2.ls))


## Selecao de modelos
AICtab(samp1.pln,samp1.ls,nobs=length(samp1),base=T, weights=T)
AICtab(samp2.pln,samp2.ls,nobs=length(samp2),base=T, weights=T)## O AICctab nao funfa, por algum problema com o nobs: VERIFICAR


##GRAFICOS##

## Amostra de uma lognormal##
## Esperados para as rank-plot com lista de coeficientes
samp1.pl.rad <- radpred(x=samp1,sad="poilog",coef=as.list(coef(samp1.pln)))
samp1.ls.rad <- radpred(x=samp1,sad="ls",coef=as.list(coef(samp1.ls)))
plot(rad(samp1))
points(samp1.pl.rad)
points(samp1.ls.rad,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

## Esperados para as oitavas com lista de coeficientes
samp1.pl.octav <- octavpred(x=samp1,sad="poilog",coef=as.list(coef(samp1.pln)))
samp1.ls.octav <- octavpred(x=samp1,sad="ls",coef=as.list(coef(samp1.ls)))
plot(octav(samp1))
points(samp1.pl.octav)
points(samp1.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

## Amostra de uma logserie##

## Esperados para as rank-plot com lista de coeficientes
samp2.pl.rad <- radpred(x=samp2,sad="poilog",coef=as.list(coef(samp2.pln)))
samp2.ls.rad <- radpred(x=samp2,sad="ls",coef=as.list(coef(samp2.ls)))
plot(rad(samp2))
points(samp2.pl.rad)
points(samp2.ls.rad,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

## Esperados para as oitavas com lista de coeficientes
samp2.pl.octav <- octavpred(x=samp2,sad="poilog",coef=as.list(coef(samp2.pln)))
samp2.ls.octav <- octavpred(x=samp2,sad="ls",coef=as.list(coef(samp2.ls)))
plot(octav(samp2),ylim=c(0,max(c(octav(samp2)$Freq,
                    samp2.pl.octav$Freq,samp2.ls.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))


## TODO: Ajuste de todos os modelos: funcao fitsad
samp2.all <- fitsad(samp2)
## TODO: graficos comparativos dos dois ajustes em cada grafico
plot(samp2.all)
