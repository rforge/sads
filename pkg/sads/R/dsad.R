dsad <- Vectorize(FUN=
                  function(y,a,lambda){
                    poi <- function(y,n){
                      w <- y*log(a*n)-lfactorial(y)-a*n
                      exp(w)
                    }
                    f1 <- function(n){
                      dexp(n,rate=lambda) * poi(y,n)
                    }
                    integrate(f1,0,Inf, abs.tol = .Machine$double.eps^0.25/1000000)$value
                  },
                  "y")
