## Tentativa de implementar o modelo de zsm proposto por Alonso & McKane 2004 Ecol Letters - Sampling Neutral ...
library(bbmle)
library(untb)
library(MASS)
library(sads)
library(vegan)
##Equacao A5: probabilidade de n individuos na amostra de tamanho J,
## dados abund proporcional de x na metacomunidade, e probabilidade de migracao m
alonsoA5 <- function(n, J, m, x){
  gama <- m*(J-1)/(1-m)
  nu <- J+gama*(1-x)
  lambda <- gama*x
  y <- lchoose(J, n)+lgamma(n+lambda)-lgamma(lambda)+lgamma(nu-n)-lgamma(nu-J)+lgamma(lambda+nu-J)-lgamma(lambda+nu)
  exp(y)
}

## Equacao 10: n esperado de especies com n individuos na amostra
alonso10 <- function(n, J, m, theta, ...){
  f1 <- function(x, N){
    alonsoA5(n=N, J, m, x)*((1-x)^(theta-1))/x
  }
  f2 <- function(N){
    integrate(f1, lower=0, upper=1, N=N,...)$value
  }
  theta*sapply(n,f2)
}

sum(alonso10(1:150, 150, 0.05, 2))

## Integrador Monte-Carlo (de Jones et al Scientific Simulation using R)
mc.integral <- function(ftn, a, b, n=1e4, ...) {
# Monte-Carlo integral of ftn over [a, b] using a sample of size n
u <- runif(n, a, b)
x <- sapply(u, ftn, ...)
return(mean(x)*(b-a))
}

## Equacao 10 com integrador Monte Carlo: n esperado de especies com n individuos na amostra
## MUITO LENTO não use
alonso10b <- function(n,J,m,theta,...){
  f1 <- function(x,N){
    alonsoA5(n=N,J,m,x)*((1-x)^(theta-1))/x
  }
  f2 <- function(N){
  mc.integral(f1,a=0, b=1,N=N,...)
}
  theta*sapply(n,f2)
}

## equacao 14: pdf
alonso14 <- function(n, J, m, theta, log = FALSE,...){
  all.values <- alonso10(1:J, J, m, theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}

alonso14 <- function(n, J, m, theta, log = TRUE,...){
  all.values <- alonso10(1:J, J, m, theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}
## Testes 
system.time((teste <- alonso10(1:100,J=100,m=0.77,theta=40)))
sum(teste)
system.time((teste2 <- alonso10b(1:100,J=100,m=0.77,theta=40)))## MUUUUITO mais lento
(alonso14(n=1:5, J=1500, m=0.77, theta=41))# problema de convergencia
plot(teste, teste2)
abline(0,1)## integrador MC tem problemas com valores altos

## Teste com um conjunto de dados
library(Hughes)
ARN82.eN.dec78

## Funcao de verossimilhanca

L1 <- function(M,T){
  -sum(alonso14(ARN82.eN.dec78,J=sum(ARN82.eN.dec78),m=M,theta=T,log=TRUE, stop.on.error=FALSE))# nao para com erros de integracao
}
lk <- logkda.R(count(ARN82.eN.dec78), use.brob=TRUE)
optimal.params(D=count(ARN82.eN.dec78), log.kda=lk)
teste <- mle2(L1,start=list(M=0.67,T=6.1)) ## problemas de convergencia, da integracao
## Tentativa com bounded values e sem hessiana
teste <- mle2(L1,start=list(M=0.67,T=6.1), skip.hessian=TRUE,method="L-BFGS-B",
              lower=c(0.001,2), upper=c(0.9999,20)) ## problemas de convergencia, da integracao
summary(teste) ## funfou!
logLik(teste)
cf <- coef(teste)

## Grafico de oitavas
library(sads)
plot(octav(ARN82.eN.dec78))
previsto <- alonso10(n=1:max(ARN82.eN.dec78),J=sum(ARN82.eN.dec78),m=cf[1],theta=cf[2])
(oitavas <- cut(c(1:length(previsto)),breaks=c(0,2^(0:10))))
plot(octav(ARN82.eN.dec78))
points(0.5:10.5,xtabs(previsto~oitavas), type="b")
## Lognormal
A.ln <- fitpoilog(ARN82.eN.dec78)
(cf2 <- as.list(coef(A.ln)))
AIC(A.ln)
AIC(teste)
points(octavpred(ARN82.eN.dec78,sad="poilog",coef=cf2), col="blue") ## ok
## Superficie de verossimilhanca
L1(M=0.67,T=6.1)
L1(M=0.6740454,T=6.1155117)
L1(M=0.5,T=2)
L1v <- Vectorize(L1)
logLik(fitdistr(ARN82.eN.dec78, "lognormal")) ## compativel!
M <- seq(0.05,0.7, by=0.025)
T <- seq(1,10,by=0.25)
m1 <- outer(M,T,"L1v")
save.image()
contour(M,T,m1-min(m1))
image(M,T,m1-min(m1))
## parece ser o ponto de menor LL
L1(M=0.4244882,T=6.758624)

## Verificando com os dados do Fisher et al 1943
data(moths)
sum(moths)
length(moths)
L2 <- function(M,T){
  -sum(alonso14(moths, J=sum(moths),m=M, theta=T, log=TRUE,
                rel.tol = sqrt(.Machine$double.eps), 
                subdivisions = 500))
}

## com bounded values e sem hessiana
teste2 <- mle2(L2, start=list(M=0.77,T=41), skip.hessian=TRUE, method="L-BFGS-B",
              lower=c(M=0.001,T=20), upper=c(M=0.9999,T=50)) 
summary(teste2) ## valor de m bem diferente do reportado por alonso

AIC(teste2)
## Comparacao com valores obtidos com o Tetame
logLik(teste2)
(teste2.LL.tet <- L2(M=0.519079,T=43.1623))
logLik(teste2)+teste2.LL.tet ## pequena diferenca
## Resultados do Tetame
teste2.tetame <- read.table("fisher43_out.txt")
(teste2.tetame)


## Ajuste aa logserie
teste2.ls <- fitlogser(moths)
summary(teste2.ls)
teste2.rad <- radpred(moths, sad="ls", coef=as.list(coef(teste2.ls)))
plot(rad(moths))
points(teste2.rad)
## Comparacao dos modelos: mto parecidos
AICtab(teste2,teste2.ls)
## Com o untb nao funfa
lkda2 <- logkda.R(count(moths), use.brob=TRUE)## resulta em algo que nao funciona com optimal params
lkda2 <- logkda(count(moths)) ## resulta em NA´ s
optimal.params(D=count(moths), log.kda=lkda2,start=c(41, 0.77) ) ## nao funfou nem com starting values
