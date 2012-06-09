There is a huge bias in the results returned by dpoith and the results returned by dsad that integrates the compound
distribution of a poisson sampling of a truncated hyperbolic.

From the single example in the logbook:

## bias expressed as differences relative to the value retuened by dpoith
library(HyperbolicDist)
library(zipfR)
y <- 0:150
y1 <- dsad(y,frac=0.05,hyperb,samp="Poisson",Theta=c(2,2,2,2))
y2 <- dpoith(y,frac=0.05)
x <- y1-y2/y2
plot(y,x,ylab="(dsad - dpoith)/dpoith")
