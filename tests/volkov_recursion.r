## Construncting a list with the terms of the recursion for solvind eq 21 in Box 2 of Volkov et al 2007
vlist <- vector(mode="list",length=10)
vlist[1] <- expression(1/omega)
cat("1/omega", file="temp.txt")
for(i in 2:10){
df <- D(vlist[[i-1]],"omega")
f.text <- scan("temp.txt",what="character",sep=";")
cat("(",f.text,")*",i-1,"-(",deparse(df),")", file="temp.txt",fill=TRUE)
vlist[i] <- parse("temp.txt")
}
## Verificando
## Novamente BCI, dados e ajustes  Volkov
J=21457;theta=48;m=0.09
bci.pred <- v2(n=1:10,theta,m,J)
m.til <- J*m/(1-m)
omega <- theta/m.til-log(m)
## n=1
v.f2(1,m,theta)*eval(vlist[[1]])
bci.pred[1]#ok
## n=3
v.f2(3,m,theta)*eval(vlist[[3]])
bci.pred[3]
## n=10
v.f2(10,m,theta)*eval(vlist[[10]])
bci.pred[10] ## funfou, mas demora muitissimo para construir a lista

## Tentando com Ryacas
library(Ryacas)
yacas(expression(1/y))
x <- Sym("x")
(f1 <- deriv(1/x,x))
for(i in 2:3){
  f2 <- deriv(f1,x)
  f3 <- (i-1)*f1-f2
  print(f3)
  f1 <- f3
}
