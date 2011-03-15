dsad_gamma <- Vectorize(FUN=
                  function(y,a,lambda,shape=0.75,fn=2){
                    poi <- function(y,n){
                      w <- y*log(a*n)-lfactorial(y)-a*n
                      exp(w)
                    }
                    if(fn==1){
                      f1 <- function(n){
                        dexp(n,rate=lambda) * poi(y,n)
                      }
                    }else{					
                      f1 <- function(n){
                        dgamma(n, shape, rate = lambda, scale = 1/rate) * poi(y,n)
                      }
                    }
                    if (y<10){
                    	integrate(f1,0,100/max(a,lambda))$value
                    }else if (y<50){
                    	k <- 6
                    	integrate(f1,0,k*y/(a+lambda))$value
                    }else{
                    	k <- 1.9920941
                    	integrate(f1,0,k*y/(a+lambda))$value
                    }
                  },
                  "y")
