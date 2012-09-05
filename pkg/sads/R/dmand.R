dmand <- function (x, N, s, v, log = FALSE) {
  if (any(x < 1)) 
    warning("the Zipf-Mandelbrot's distribution is not set to zero or negative numbers")
  if (s <= 0)
    stop("s must be greater than zero")
  if (v < 0) 
    stop("v must be greater than zero")
  if (N < 1) 
    stop("N must be positive integer")
  if (!any(is.wholenumber(x))) 
    warning("x must be integer")
  lny <- - s * log(x+v) - log(sum(1/((1:N)+v)^s))
  if (log) return(lny)
  else return(exp(lny))
}