##There is a huge bias in the  results returned by the dpoig function in comparison with the results returned by dsad with the same parameters.

##See the single example taken from the logbook below, biases expressed as the difference among values relative to the value retirned by dpoig:

## Code modified from the logbook
y <- 0:150
y1 <- dsad(y,frac=0.05,gamma,samp="Poisson",rate=1/1000,shape=1)
y2 <- dpoig(y, frac=0.05, rate=1/1000, shape=1)
x <- (y1-y2)/y2
plot(y,x,ylab="(dsad - dpoig)/dpoig")
