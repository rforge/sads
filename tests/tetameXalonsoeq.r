## Tentativa de implementar o modelo de zsm proposto por Alonso & McKane 2004 Ecol Letters - Sampling Neutral ...
## Log-likelihoods for the ZSM from Tetame and from Alonso & McKane 2004 (Ecol Lett 7: 901-910)
library(vegan)
library(bbmle)
library(MASS)
library(poilog)

## A&M equation A4:
## Probability of a given abundance n in a sample of size J for a species that
## has proportional abundance x in a metacommunity and migration probability m
alonsoA4 <- function(n,J,m,x){
  gama <- m*(J-1)/(1-m)
  nu <- J+gama*(1-x)
  lambda <- gama*x
  y <- lchoose(J,n)+lgamma(n+gama*x)-lgamma(gama*x)+lgamma(nu-n)-lgamma(nu-J)+lgamma(lambda+nu-J)-lgamma(lambda+nu)
  exp(y)
}

## A&M equation 10: expected number of species with abundance =n in a sample of size J
## This the product of equation A4 and eqution 9, integrated from zero to 1
alonso10 <- function(n,J,m,theta,...){
  f1 <- function(x,N){
    alonsoA4(n=N,J,m,x)*((1-x)^(theta-1))/x
  }
  f2 <- function(N){
  integrate(f1,lower=0, upper=1, N=N,...)$value
}
  theta*sapply(n,f2)
}
## equation 14: pdf for a given abundance under the ZSM
alonso14 <- function(n,J,m,theta, log=FALSE,...){
  all.values <- alonso10(1:J,J,m,theta,...)
  lprobs <- log(all.values[n])-log(sum(all.values))
  if(log) lprobs
  else exp(lprobs)
}

## And now the fitting by maximun likelihood, using bbmle package
## Fisher et al 1943 moth data
moths <- read.table("fisher43.txt")[[1]]
sum(moths)
length(moths)
## Log-likelihood function
## two free parameters: M=migration , T = theta
L2 <- function(M,T){
  -sum(alonso14(moths,J=sum(moths),m=M,theta=T,log=TRUE,
                rel.tol = sqrt(.Machine$double.eps), 
                subdivisions = 500))
}

## Fitting with bounded optmization to keep M between zero and one
## Takes a little time
teste2 <- mle2(L2,start=list(M=0.77,T=41), skip.hessian=TRUE,method="L-BFGS-B",
              lower=c(M=0.001,T=20), upper=c(M=0.9999,T=50))
##Fitted model
## Note that estimate of M is quite different from those reported by A&M: 0.59 here bu 0.77 in A&M
summary(teste2)

## Results from Tetame: very similar estimates, (M=0.52 and theta=43.2)
teste2.tetame <- read.table("fisher43_out.txt", header=T)
teste2.tetame
## In fact the loglikelihood above using estimates for Tetame is very close to the optmized log-likelihood
logLik(teste2)
(teste2.LL.tet <- L2(M=teste2.tetame$m,T=teste2.tetame$Theta))
## But the log-likelihoods differ
logLik(teste2)
teste2.tetame$loglike_min

## The same data fitted to a lognormal
teste3 <- poilogMLE(moths)
teste3$logLval ## same magnitude of logLik calculated by A&M formulas

## And to lognormal
## With fitdistr
logLik(fitdistr(moths, "lognormal"))
## With vegan
logLik(rad.lognormal(moths)) ## different from both log-likelihood values
