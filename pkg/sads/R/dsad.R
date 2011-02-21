dsad <- Vectorize(FUN=
                  function(y,a,lambda){
                    poi <- function(y,n){
                      w <- y*log(a*n)-lfactorial(y)-a*n
                      exp(w)
                    }
                    f1 <- function(n){
                      dexp(n,rate=lambda) * poi(y,n)
                    }
                    t1 <- integrate(f1,0,Inf)$value
		    y <- 0
                    t0 <- integrate(f1,0,Inf)$value
		    return(t1/(1-t0))
                  },
                  "y")
