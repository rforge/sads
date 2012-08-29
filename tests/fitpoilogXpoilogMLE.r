library(sads)
library(untb)


## Semente de n aleatorios
set.seed(42)# coloque qq valor, apenas para q eu possa reproduzir
## Uma mostra Poisson de uma lognormal, zeros excluídos
set.seed(42)
samp1 <- rsad(100,frac=0.2,sad=lnorm,samp="Poisson",meanlog=3,sdlog=2) # veja help da funcao rsad
## Ajuste a Poilog pela fitpolig
fit1 <- fitpoilog(samp1, trunc=0)
## ajuste pelo poilogMLE
fit2 <- poilogMLE(samp1)
## Comparacao das verosimilhancas: diferentes
logLik(fit1)
fit2$logLval
## Coeficientes: diferentes também
coef(fit1)
fit2$par
