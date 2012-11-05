qmand <- function (p, N, s, v, lower.tail = TRUE, log.p = FALSE) {
  if (s <= 0) 
    stop("s must be greater than zero")
  if (v <= 0) 
    stop("v must be greater than zero")
  if (N < 1||!any(is.wholenumber(N))) 
    stop("N must be positive integer")
  d <- NULL
  if (log.p) 
    p <- exp(p)
  if (!lower.tail) 
    p <- 1 - p
  busca <- function(U1, U2) {
    repeat {
      tt <- pmand(U2, N, s, v)
      U2 <- ifelse(tt >= U1, U2 - 1, U2)
      a1 <- U2
      U2 <- ifelse(tt < U1, U2 + 1, U2)
      a2 <- U2
      if (pmand(min(a1, a2), N, s, v) < U1 & U1 <= pmand(max(a1, a2), N, s, v)) {
        return(max(a1, a2))
      }
    }
  }
  for (i in 1:length(p)) {
    U1 <- p[i]
    U2 <- round(runif(1, min = 1, max = N))
    if (U1 <= pmand(1, N, s, v)) {
      d[i] <- 1
    }
    else if (U1 >= 1) {
      d[i] <- N
    }
    else {
      d[i] <- busca(U1, U2)
    }
  }
  return(d)
}
