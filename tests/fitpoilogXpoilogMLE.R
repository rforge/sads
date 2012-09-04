## Comparison between fitpoilog fo sads package and poilogMLE form poilog package
## Poisson sample from a lognormal (zeroes ommited)
set.seed(1913)
samp1 <- rsad(100, frac=0.15, sad=lnorm, samp="Poisson", meanlog=3, sdlog=2)

## model fit
samp1.pln1 <- fitpoilog(samp1, trunc=0)
samp1.pln2 <- poilogMLE(samp1)
## Log-likelihoods
logLik(samp1.pln1)
samp1.pln2$logLval
## Coeficients
coef(samp1.pln1)
samp1.pln2$par
