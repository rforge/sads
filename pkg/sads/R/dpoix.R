dpoix <- function(x, a, lambda, log=FALSE) {
	  b <- a^(x-1)
	  m <- lambda*b
	  n <- (lambda + a)^x
          if(log)log(m/n) else m/n
        }
