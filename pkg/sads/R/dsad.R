dsad <- Vectorize(FUN=
                  function(y,a,lambda){
                    poi <- function(y,n){
                      w <- y*log(a*n)-lfactorial(y)-a*n
                      exp(w)
                    }
                    f1 <- function(n){
                      dexp(n,rate=lambda) * poi(y,n)
                    }
					k <- 1.9920941
                    integrate(f1,0,k*y/(a+lambda))$value
                  },
                  "y")
		  
