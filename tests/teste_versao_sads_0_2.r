## Para instalar o pacote sads (se o arquivo tar.gz estiver em outro diretorio indique no primeiro argumento)
install.packages("sads_0.2.03.tar.gz", repos=NULL)# funcionou em linux, creio q nao em outros sistemas
## Carregue os pacotes
library(sads)
library(untb)

## Ajustes e graficos passo a passo, com amostras simuladas de duas comunidades

## 1. Geramos uma amostra Poisson de duas comunidades: uma lognormal e uma logserie
## Semente de n aleatorios
set.seed(42)# coloque qq valor, apenas para q eu possa reproduzir
## Uma mostra Poisson de uma lognormal
samp1 <- rsad(100,frac=0.2,sad=lnorm,samp="Poisson",meanlog=3,sdlog=2) # veja help da funcao rsad
## O mesmo com uma amostra da logserie com mesma riqueza e total de individuos
samp2 <- fisher.ecosystem(N=sum(samp1),S=length(samp1),nmax=sum(samp1))

## 2. Ajustes dos modelos de sad a cada amostra
## Amostra 1
## Ajuste a Poilog
samp1.pln <- fitpoilog(samp1)
cf.pln <- as.numeric(coef(samp1.pln))# extracao dos coeficientes, para uso abaixo
summary(samp1.pln) # resumo do modelo, ha tb metodos para intervalos, etc, veja help da mle2
## Ajuste a Logserie
samp1.ls <- fitlogser(samp1)
cf.ls <- as.numeric(coef(samp1.ls))
## Amostra 2
## Poilog
samp2.pln <- fitpoilog(samp2)
cf2.pln <- as.numeric(coef(samp2.pln))
## ajuste à logserie
samp2.ls <- fitlogser(samp2)
cf2.ls <- as.numeric(coef(samp2.ls))

## 3. Comparamos os modelos com AIC
AICtab(samp1.pln,samp1.ls,nobs=length(samp1),base=T, weights=T)# Poilog deve vencer
AICtab(samp2.pln,samp2.ls,nobs=length(samp2),base=T, weights=T) # Logserie deve vencer

## 4. Esperados para graficos de rank-abundancia

## Para a amostra de uma lognormal
## Argumentos da funcao radpred: x=vetor de abundancias, sad=modelo de sad ("poilog" ou "ls"), coef= lista dos coeficientes do modelo de sad ajustado com poilog fit ou fitlogser
## Previsto pela poilog
samp1.pl.rad <- radpred(x=samp1,sad="poilog",coef=as.list(coef(samp1.pln))) # pode demorar
## O resultado da funcao radpred eh um dataframe com os ranks e os previstos
samp1.pl.rad
## Previsto pela logserie
samp1.ls.rad <- radpred(x=samp1,sad="ls",coef=as.list(coef(samp1.ls)))
## Graficos de observados e previstos calculados acima
plot(rad(samp1)) # a funcao rad gera um dataframe com ranks e abundancias observados. O plot aplicado sobre o objeto resultante ja entende que eh para fazer um radplot
points(samp1.pl.rad)
points(samp1.ls.rad,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))

## Amostra de uma logserie
## Esperados para as rank-plot com lista de coeficientes
samp2.pl.rad <- radpred(x=samp2,sad="poilog",coef=as.list(coef(samp2.pln))) # pode demorar
samp2.ls.rad <- radpred(x=samp2,sad="ls",coef=as.list(coef(samp2.ls)))
plot(rad(samp2))
points(samp2.pl.rad)
points(samp2.ls.rad,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))


## 5. Esperados para grafico de oitavas
## Esperados para as oitavas com lista de coeficientes
samp1.pl.octav <- octavpred(x=samp1,sad="poilog",coef=as.list(coef(samp1.pln)))
samp1.ls.octav <- octavpred(x=samp1,sad="ls",coef=as.list(coef(samp1.ls)))
plot(octav(samp1))
points(samp1.pl.octav)
points(samp1.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))


## Esperados para as oitavas com lista de coeficientes
samp2.pl.octav <- octavpred(x=samp2,sad="poilog",coef=as.list(coef(samp2.pln)))
samp2.ls.octav <- octavpred(x=samp2,sad="ls",coef=as.list(coef(samp2.ls)))
plot(octav(samp2),ylim=c(0,max(c(octav(samp2)$Freq,
                    samp2.pl.octav$Freq,samp2.ls.octav$Freq))))
points(samp2.pl.octav)
points(samp2.ls.octav,col="red")
legend("topright", c("Poisson-lognormal", "Logseries"), lty=1, col=c("blue","red"))
