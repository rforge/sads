## Teste da funcao trueLL
library(sads)
library(MASS)
## Passo a passo:
x <- rpois(10,10)
xu <- plnorm(x+0.5,meanlog=mean(log(x)),sdlog=sd(log(x)))
xl <- plnorm(x-0.5,meanlog=mean(log(x)),sdlog=sd(log(x)))
sum(log(xu-xl))
trueLL(x,lnorm,meanlog=mean(log(x)),sdlog=sd(log(x))) # ok!
##
x <- round(rlnorm(500,meanlog=3,sdlog=0.5),1)
logLik(fitdistr(x,"lognormal"))
sum(dlnorm(x,mean(log(x)),sd(log(x)),log=T))
trueLL(x,lnorm,precision=0.1,meanlog=mean(log(x)),sdlog=sd(log(x)))
## Para precisao = 1
x <- round(rnorm(500,3,5),0)
logLik(fitdistr(x,"normal"))
sum(dnorm(x,mean(x),sd(x),log=T))
trueLL(x,norm,precision=1,mean=mean(x),sd=sd(x))
sum(dnorm(x,0,1,log=T))
trueLL(x,norm,precision=1,mean=0,sd=1)
