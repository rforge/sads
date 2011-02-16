depois <- function(x, a=0.1, lambda=0.05) {

	  poix <- function(n) {

		    b <- lambda*a^(x)/factorial(x)
		    c <- -n*(a+lambda)
		    m <- exp(c)*n^(x)
		    b*m
	  }
	  begin <- 0
	  end <- Inf
	  t1 <- integrate(poix, begin, end)
	  x <- 0
	  t0 <- integrate(poix, begin, end)
	  t <- t1$value/(1-t0$value)
	  return(t)
}
