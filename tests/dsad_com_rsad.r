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

## Um teste qquick and dirty
dsadr <- function(x,N=1e5,...){
  y <- rsad(N,...)
  z <- length(y)/N
  new.N <- (N-length(y))/z
  y <- c(y,rsad(new.N,...))
  ty <- table(y)
  Fn <- approxfun(x=as.numeric(names(ty)), y=ty/sum(ty), method="constant", f=1, rule=2)
  py <- Fn(x)
  py[py==0] <- 1/sum(ty) # gambiarra para nao dar probabilidades =0
}
## Teste
y <- dsadr(1:20,N=1e6,frac=0.1,sad=lnorm,samp="Poisson",meanlog=3,sdlog=1.5)
x <- dpoilog(1:20,mu=3+log(0.1),sig=1.5)
plot(x,y, log="xy"); abline(0,1) ## mto bom
## Teste de ajuste
dados <- rsad(100,frac=0.1,sad=lnorm,samp="Poisson",meanlog=3,sdlog=2)
LL <- function(m,s){
  -sum(log(dsadr(dados,frac=0.1,sad=lnorm,samp="Poisson",meanlog=m,sdlog=s)))
}
teste <- mle2(LL,start=list(m=1,s=1))
teste2 <- poilogMLE(dados)
coef(teste)
teste2$par
## Nao funfou: demora, dá erros e estimativas muito fora
