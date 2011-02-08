dpoix <- function(x, a = 0.1, lambda = 0.5) {
	  
	  b <- a^(x-1)
	  m <- lambda*b
	  n <- (lambda + a)^x
	  m/n
}