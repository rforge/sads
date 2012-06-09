samp1 <- rsad(100,frac=0.1,sad=lnorm,samp="Poisson",meanlog=3,sdlog=2)

samp1.fit <- fitpoilog(samp1)
cf <- as.numeric(coef(samp1.fit))
t1 <- radpred(samp1,"poilog",mu=cf[1],sig=cf[2])
t2 <- radpred(samp1,"poilog",mu=3+log(0.25),sig=2)
radplot(samp1)
lines(t1)
lines(t2,col="blue")

samp2 <- fisher.ecosystem(N=sum(samp1),S=length(samp1),nmax=sum(samp1))
t3 <- radpred(samp2,"ls",N=sum(samp2),
              alpha=fishers.alpha(N=sum(samp2),S=length(samp2)))
radplot(samp2)
lines(t3)

histsad(samp1,prob=T) ## apague esta funcao
t1b <- hpred(samp1,"poilog",mu=cf[1],sig=cf[2]) ## e esta tb
points(t1b,type="b")

## Funcao com metodo S4
cf <- as.numeric(coef(samp1.fit))
samp1.oct <- octav(samp1)
plot(samp1.oct)
octavpred(x=samp1,sad="poilog",mu=cf[1],sig=cf[2])
points(octavpred(x=samp1,sad="poilog",mu=cf[1],sig=cf[2]))
