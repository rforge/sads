dsad <- Vectorize(FUN=
                  function(y,a,lambda){
                    poi <- function(y,n){
                      w <- y*log(a*n)-lfactorial(y)-a*n
                      exp(w)
                    }
                    f1 <- function(n){
                      dexp(n,rate=lambda) * poi(y,n)
                    }
					inf <- 5000
                    integrate(f1,0,inf, abs.tol = .Machine$double.eps)$value
                  },
                  "y")
