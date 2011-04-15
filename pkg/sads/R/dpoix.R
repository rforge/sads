dpoix <- function(y, frac, rate, log=FALSE) {  
	  b <- y*log(frac)
	  m <- log(rate)
	  n <- (y+1)*log(rate+frac)
    if(log)b+m-n else exp(b+m-n)
}
y <- 0:50
y1 <- dsad(y,frac=0.05,gamma,samp="Poisson",rate=1/1000,shape=2)
y2 <- dpoith(y=1,frac=0.05)
x <- ((y1-y2)^2)^.5
plot(y,x,ylab="|dsad|-|dpoig|")