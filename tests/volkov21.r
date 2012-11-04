## Tentativa de implementar a funcao 21 de Volkov
v2 <- function(n, theta, m, J){
  X <- 1-m
  m.tilde <- J*m/(1-m)
  omega <- theta/m.tilde - log(1-X)
  abunds <- unique(n)
  f1 <- function(y,N){
    z <- (lgamma(N+y)-lgamma(1+y))-omega*y
    exp(z)
  }
  f2 <- function(ab){
    k <- exp(log(theta)+ab*log(X)-lfactorial(ab))
    k*integrate(f1,0,Inf, N=ab)$value
  }
  vals <- sapply(abunds,f2)
  names(vals) <- abunds
  vals[as.character(n)]
}

## Previstos para primeiras oitavas em BCI, usando mles da tabela 1
(bci.v <- v2(n=1:8,theta=48,m=0.09,J=21457))
## parece que a figura esta usando o metodo do Preston de redistribuir metade dos valores para cada lado
bci.v2 <- bci.v
bci.v2[1] <- bci.v[1]/2
bci.v2[2] <- bci.v[2]+bci.v[1]/2-bci.v[2]/2
bci.v2 ## parece que e isto
## Oitavas sem essa maracutaia
data(BCI)
bci <- apply(BCI,2,sum)
plot(octav(bci))
(bci.v <- v2(n=1:8,theta=48,m=0.09,J=21457))
## terceira oitava: spuerestima um pouco
sum(bci.v[3:4])
## quarta oitava: superestima, mas compativel
sum(bci.v[5:8])

################################################################################
## A integral nao funciona para valores muito grandes
## Mas VOlkov indica no Box 2 que a integral f(n,omega) pode ser reolvida por recursao:
## f(1,omega)=1/omega
## f(n+1,omega)=n*f(n,omega)-D(f(n.omega),"omega")

## Conferindo para os dois primeiros valores
## Novamente BCI, dados e ajustes  Volkov
J=21457;theta=48;m=0.09
bci.pred <- v2(n=1:3,theta,m,J)
m.til <- J*m/(1-m)
omega <- theta/m.til-log(m)
## Funcao para o termo que multiplica integral na expressao 21
v.f2 <- function(n,m,theta){
  exp(log(theta)+(n*log((1-m))-lfactorial(n)))
}
##Para n=1 o valor previsto sera
v.f2(1,m,theta)*1/omega ## ok!
## para o segundo temos
f1 <- expression(1/omega)
df1 <- D(f1,"omega")
temp <- paste("1/omega-",deparse(df1))
f2 <- parse(text=temp)
v.f2(2,m,theta)*eval(f2) ## correto
## Para n=3
df2 <- D(f2,"omega")
temp <- paste("(",temp,")*2-",deparse(df2))
f3 <- parse(text=temp)
v.f2(3,m,theta)*eval(f3) ## correto
## Rapidamente esta iteracoes resultam em expressoes enormes!

################################################################################
## Estudo da funcao que e a integral
f21a <- function(y,omega,n){
  z <- (lgamma(n+y)-lgamma(1+y))-omega*y
  exp(z)
}

f21b <- function(theta,m,n){
  z <- log(theta)+(n*log(1-m)-lfactorial(n))
  exp(z)
}

f21b(theta,m,100)

## Uma funcao com tudo junto
f21c <- function(y,omega,theta,m,n){
  z <-  log(theta)+(n*log(1-m)-lfactorial(n)) + ((lgamma(n+y)-lgamma(1+y))-omega*y)
  exp(z)
}
f21c(1000,omega,theta,m,1000)
integrate(f21c,lower=0,upper=Inf,omega=omega,theta=theta,m=m,n=5000)
## Para um individuo
integrate(f21c,lower=0,upper=Inf,omega=omega,theta=theta,m=m,n=1)
bci.pred
## Para dois
integrate(f21c,lower=0,upper=Inf,omega=omega,theta=theta,m=m,n=2)## ok!


## Nova versao da funcao, com constante para dentro da integral
v2b <- function(n, theta, m, J, log=FALSE, tol=1e-4, ...){
  X <- 1-m
  m.tilde <- J*m/(1-m)
  omega <- theta/m.tilde - log(1-X)
  #abunds <- 1:(J/10)
  f1 <- function(y,N){
    k <- log(theta)+(N*log(X)-lfactorial(N)) ## parte constante q esta fora da integral na eq 21
    f <- (lgamma(N+y)-lgamma(1+y))-omega*y ## parte q esta do lado de dentro
    exp(k+f)
  }
  f2 <- function(ab){
   integrate(f1,0,Inf, N=ab, ...)$value
  }
  abunds <- 1
  vals <- c()
  vals[1] <- f2(abunds)
  abunds <- abunds+1
  vals[2] <- f2(abunds)
  while((sum(vals))/(sum(vals[-abunds]))>(1+tol)|abunds<max(n))
    {
      abunds <- abunds+1
      vals[abunds] <- f2(abunds)
      }
  ES <- theta*log((1-J*m)/(theta*(1-m))*log(m))
  if(ES>sum(vals)) Stot <- ES
  else Stot <- sum(vals)
  if(log)log(vals[n]/Stot) else vals[n]/Stot
}

v2b(max(bci), theta=theta,m=m,J=J, tol=1e-6, log=T)
## Verossimilhanca de bci
-sum(v2b(bci,theta,m,J=sum(bci),tol=1e-4,log=T))
## Comparando com poilog
bci.pl <- fitsad(bci,"poilog")
logLik(bci.pl) ## valores na mesma escala
## Com moths, usando estimativas de Alonso & McKane
moths.zsm <- optimal.params(moths, moths.kda) ## com estes nao funfa
moths.zsm
-sum(v2b(moths,theta=41,m=0.77,J=sum(moths),tol=1e-4,log=T)) ## dah pau mesmo com estimativas do Alonso
## sem o maior valor
temp <- moths[!moths==max(moths)]
-sum(v2b(temp,theta=41,m=0.77,J=sum(temp),tol=1e-4,log=T)) ## agora foi, mas loglik altissima
## Tantando o untb
temp.kda <- logkda.R(temp)
temp.zsm <- optimal.params(temp, temp.kda) ## com estes nao funfa
temp.zsm ## m muito alto
-sum(v2b(temp,theta=temp.zsm[1],m=temp.zsm[2],J=sum(temp),tol=1e-4,log=T)) ## deu pau
## Logseries
temp.ls <- fitsad(temp,"ls")
logLik(temp.ls)
## Funcao para encontrar valores iniciais baseada nas eqs 6 e 7
fv6 <- function(x,S,J){
  #theta <- exp(x[1])
  #m <- exp(x[2])/(1+exp(x[2]))
  theta <- x[1]
  m <- x[2]
  abs(S-theta*log(1-J*m/(theta*(1-m))*log(m)))
}
teste <- optim(c(45,0.5),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(c(4,0.9),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(c(10,0.9),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
## Fazendo iteracoes
teste <- optim(c(10,0.5),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par ## MUITO instavel
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="SANN")
teste$par
## partindo do alfa de fisher
st <- coef(fitsad(bci,"ls"))
teste <- optim(c(st,0.99),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
## partindo od theta se m=1
st <- optimal.theta(bci)
teste <- optim(c(st,0.5),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par 
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par ## MUITO instavel
teste <- optim(jitter(teste$par),fv6,S=225,J=21457, method="SANN")
teste$par
## A mesma funcao com o m fixo
fv6b <- function(x,S,J,m){
  theta <- x
  S-theta*log(1-J*m/(theta*(1-m))*log(m))
}
## Com uniroot para calcular apenas theta nos extremos de m e depois usando a media desses valores com m=0.5
cf1 <- uniroot(fv6b,c(0.1,5000),S=225,J=21457,m=0.01)$root
cf2 <- uniroot(fv6b,c(0.1,5000),S=225,J=21457,m=0.99)$root
teste <- optim(c(mean(c(cf1,cf2)),0.5),fv6,S=225,J=21457, method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
## Boa! Tentando com moths
(cf1 <- uniroot(fv6b,c(0.1,5000),S=length(moths),J=sum(moths),m=0.01)$root)
(cf2 <- uniroot(fv6b,c(0.1,5000),S=length(moths),J=sum(moths),m=0.999)$root)
teste <- optim(c(mean(c(cf1,cf2)),0.9),fv6,S=length(moths),J=sum(moths), method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par
## Verificando com untb
moths.kda <- logkda.R(moths)
(cf3 <- optimal.params(moths, log.kda=moths.kda)) ## nao muito bom
## usando estes como valor de partida
teste <- optim(cf3,fv6,S=length(moths),J=sum(moths), method="L-BFGS-B",lower=c(0.1,0.000001),upper=c(1000,0.999999))
teste$par ## mesmo assim foge muito
coef(fitsad(moths, "ls"))
## Com funcoes de ligacao (piora)
fv6b <- function(x,S,J){
  theta <- exp(x[1])
  m <- exp(x[2])/(1+exp(x[2]))
  abs(S-theta*log(1-J*m/(theta*(1-m))*log(m)))
}
teste <- optim(c(log(st),log(0.5/0.5)),fv6b,S=225,J=21457, method="SANN")
exp(teste$par[1])
exp(teste$par[2])/(1+exp(teste$par[2]))
