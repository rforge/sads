library(untb)
library(bbmle)
library(vegan)

## Densidade a partir da funcao do Volkov
dvolkov <- function(x, J, m, theta, log=T)
{
  Sj <- volkov(J=J,params=c(theta=theta,m=m))
  y <- Sj[x]/sum(Sj)
  if(log)y else y
}

dvolko(1:5,J=21386,theta=50,m=0.1)
teste3 <- dvolko(1:21386,J=21386,theta=10,m=0.00001)
teste2 <- volkov(J=21386,params=c(theta=10, m=0.00001))
teste3[1:5]
p2[1:5]
p2 <- teste2/sum(teste2)
-sum(log(p2[bci])) ## logverossimilhanca
## igua a
LL(theta=10,m=0.00001)
LL(theta=44,m=0.15)## valores do McKane
LL(theta=95.8,m=0.002)## valores com a funcao dzsm antiga
## Comparando com dvolko
teste2b <- dvolko(1:sum(bci),J=sum(bci),theta=10,m=0.00001)
teste2b[1:5]
teste <- volkov(J=21386,params=c(theta=50, m=0.1))
teste[1:5]/sum(teste)

## Teste com BCI
data(BCI)
bci <- apply(BCI,2,sum)
LL <- function(theta,m)
  {
    -sum(log(dvolko(bci,J=sum(bci),theta=theta, m=m)))
  }
teste.bci <- mle2(LL, start=list(theta=50, m=0.1), method="L-BFGS-B", lower=c(40,0.001), upper=c(100,0.9))## problema
summary(teste.bci)


## Com logkda e optimal.params do untb
moths.kda <- logkda.R(moths)
moths.zsm <- optimal.params(moths, log.kda=moths.kda)
## Testes da nova funcao de zsm pela formula do volkov, usando parametros obitidos acima
moths.volk <- fitvolkov(moths, start.value=c(0.99, 50))

