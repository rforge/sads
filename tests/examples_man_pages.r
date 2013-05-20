## fitbs ##
## Magurran (1989) example 5:
## birds in an Australian forest
mag5 <- c(103,115,13,2,67,36,51,8,6,61,10,21,
          7,65,4,49,92,37,16,6,23,9,2,6,5,4,
          1,3,1,9,2)
mag5.bs <- fitbs(mag5)
summary(mag5.bs)## no estimated coefficient
coef(mag5.bs) ## fixed coeficients N and S
## Diagnostic plots
par(mfrow=c(2,2))
plot(mag5.bs)
par(mfrow=c(1,1))

### mzsm ###
## Alonso & McKanne (2004) fig 2
data(moths) #Fisher's moths data
m.tab <- hist(moths, breaks = 2^(0:12), plot = FALSE)
plot(m.tab$density~m.tab$mids, log="xy",
     xlab = "Abundance", ylab = "Probability density",
     ylim=c(1e-7,1))
X <- 1:max(moths)
Y <- dmzsm(X, J = sum(moths), theta = 39.8)
lines(Y ~ X)

### fitmzsm ###
data(moths) #Fisher's moths data
moths.mzsm <- fitmzsm(moths) ## same as fitsad(moths, sad="mzsm")
coef(moths.mzsm) ## Compare with theta=38.9, Alonso&McKanne (2004)
## Diagnostic plots
par(mfrow=c(2,2))
plot(moths.mzsm)
par(mfrow=c(1,1))
## Comparison to logseries
## fit to logseries
moths.ls <- fitsad(moths, sad="ls")
## Model selection
AICtab(moths.ls, moths.mzsm, weights=T)

## dtrunc, ptrunc ##
A <- dtrunc("lnorm", x = 1:5, trunc = 0.5,
       coef = list( meanlog=1, sdlog=1.5 ) )
## same as
B <- dlnorm( 1:5 , meanlog = 1, sdlog = 1.5 ) /
  ( plnorm ( 0.5 , meanlog = 1, sdlog = 1.5, lower = FALSE))
## checking
identical( A, B )

A <- ptrunc("pois", q = 1:5, trunc = 0,
       coef = list( lambda = 1.5 ) )
## same as
B <- (ppois( 1:5 , lambda = 1.5 ) -
      ppois(0 , lambda = 1.5 ) ) /
  (ppois(0 , lambda = 1.5, lower = FALSE))
## checking
identical(A,B)


## dvolkov ##
## Volkov et al 2003 fig 1
## (only fit ro Volkov's model
library( vegan )
data( BCI )
bci.oct <- octav( bci, preston=TRUE )
plot( bci.oct )
cdf <- pvolkov( bci.oct$upper, theta = 47.226, m = 0.1, J = sum(bci) )
bci.exp <- diff( c(0,cdf) ) * length(bci)
midpoints <- as.numeric( bci.oct$octave ) - 0.5
lines( midpoints, bci.exp, type="b" )
## the same with octavpred
bci.exp2 <- octavpred( bci, sad = "volkov",
                      coef = list(theta = 47.226, m = 0.1) )
lines( bci.exp2 )

### fitvolkov ###
## Magurran (1989) example 5:
## birds in an Australian forest
mag5 <- c(103,115,13,2,67,36,51,8,6,61,10,21,
          7,65,4,49,92,37,16,6,23,9,2,6,5,4,
          1,3,1,9,2)
mg5.v <- fitvolkov(mag5)
coef(mg5.v)
par(mfrow=c(2,2))
plot(mg5.v)
par(mfrow=c(1,1))

## fitsad ##
moths.ln <- fitsad(moths, "lnorm", trunc=0.5)
moths.pln <- fitsad(moths, "poilog")
AICctab(moths.ln, moths.pln, nobs=length(moths), weights=TRUE)
plot(octav(moths))
lines(octavpred(moths.ln))
lines(octavpred(moths.pln), col="red")
legend()
## Biomass as abundance variable
data(ARN82.eB.apr77) #bentonic marine animals
AR.ln <- fitsad(ARN82.eB.apr77, sad="lnorm", dec.places=2)
AR.g <- fitsad(ARN82.eB.apr77, sad="gamma", dec.places=2)
AR.wb <- fitsad(ARN82.eB.apr77, sad="weibull", dec.places=2)                
plot(octav(ARN82.eB.apr77))
lines(octavpred(AR.ln))
lines(octavpred(AR.g), col="red")
lines(octavpred(AR.wb), col="green")
legend("topright", c("lognormal", "gamma", "weibull"), lty=1, col=c("blue","red", "green"))
AICctab(AR.ln, AR.g, AR.wb, nobs=length(ARN82.eB.apr77), weights=T)
##alternativa: rads
plot(rad(ARN82.eB.apr77))
lines(radpred(AR.ln))
lines(radpred(AR.g), col="red")
lines(radpred(AR.wb), col="green")

## Octavpred ##
## effect of Preston criteria
data(moths) ## Fisher's moth data
moths.ln <- fitsad(moths, "lnorm") ## fit to the lognormal
AICctab(moths.ln, moths.pln, nobs=length(moths), weights=TRUE)
par(mfrow=c(1,2))
plot(octav(moths))
lines(octavpred(moths.ln))
plot(octav(moths, preston=T))
lines(octavpred(moths.ln))
par(mfrow=c(1,1))
