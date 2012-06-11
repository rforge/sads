## Uma funcao generica de densidade a partir od sorteio de valores pela rsad

## Uma amostra de uma lognormal
am1 <- rsad(100,frac=0.1,sad=lnorm,samp="Poisson",meanlog=3,sdlog=1.5)
## Uma amostra enorme desta mesma distribuicao
set.seed(42)
d.am1 <- rsad(1e3,frac=0.1,sad=lnorm,samp="Poisson",meanlog=3,sdlog=1.5)
length(d.am1)
table(d.am1)
Fn <- ecdf(d.am1) ## cria uma funcao cdf empirica
## Ha um buraco entre 29 e 30
Fn(28)
Fn(29)
Fn(30)
## A ecdf retorna as proporcoes observadas apenas, tentando com approxfun
x <- table(d.am1)
Fn2 <- approxfun(x=as.numeric(names(x)), y=x/sum(x), method="constant", f=1) ## f 1 da uma pdf mais conservadora para likelihood, supondo uma faltarao valores na cauda direita e que esta cauda decai monotonicamente.
Fn2(28)
Fn(28)-Fn(27)
Fn(30)-Fn(29)
Fn2(1)
Fn(1)
## Comparando resultados
Fn2(28)
Fn2(29)
Fn2(30)
Fn2(29)==Fn2(30)
