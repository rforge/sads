dpoix <- function(x, frac, rate, log=FALSE) {
	  b <- frac^(x-1)
	  m <- rate*b
	  n <- (rate + frac)^x
          if(log)log(m/n) else m/n
        }
