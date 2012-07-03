## Teste da primeira versao de radpred para funcoes continuas
## Para funcionar carregue todas as funcoes do pacote sads e os pacotes dependencias
library(sads)
## A funcao
radpredc <- function(x, sad, coef,...){
  if(missing(coef))dots <- list(...)
  else dots <- c(coef, list(...))
  qsad <- get(paste("q", deparse(substitute(sad)), sep=""), mode = "function")
  S <- length(x)
  Y <- ppoints(S)## Gerando os dados
  ab <- do.call(qsad,c(list(p=Y,lower.tail=F),dots))
  new("rad",data.frame(rank=1:S, abund=ab))
}

## Dados de uma lognormal ##
## Gerando os dados
(x1 <- rlnorm(100,meanlog=0.5,sdlog=1))
##Ajuste
x1.fit <- fitdistr(x1,"lognormal")
## Pegando coeficientes
x1.cf <- coef(x1.fit)
(x1.radpred <- radpredc(x1,lnorm,coef=as.list(x1.cf)))
## Grafico rank-abund com observados e esperados
x1.rad <- rad(x1) ## gera tabela com ranks e abund observadas
plot(x1.rad) ## plota
points(x1.radpred) ## Adiciona a linha

## Segundo teste: uma weibull
## Gerando os dados
(x2 <- rweibull(100,shape=0.8,scale=3))
##Ajuste
x2.fit <- fitdistr(x2,"weibull")
## Pegando coeficientes
x2.cf <- coef(x2.fit)
(x2.radpred <- radpredc(x2,weibull,coef=as.list(x2.cf)))
## Grafico rank-abund com observados e esperados
x2.rad <- rad(x2) ## gera tabela com ranks e abund observadas
plot(x2.rad) ## plota
points(x2.radpred) ## Adiciona a linha
