dls <- function(x, N, alpha, log=F){
		X<-N/(N+alpha)
		gama<- function(y) {1/log(1/(1-y))}
		y <- gama(X)*(X^x)/x
                if(log)log(y)
                else y
              }
