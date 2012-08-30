fitlogser <- function(x, size, rich, trunc, start.value, upper = 20, ...){
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (!missing(x)){
    S <- length(x)
    N <- sum(x)
    if (!missing(size)| !missing(rich)){
      warning(paste("Model fitted with size = ", N, " and rich = ", S, " \n calculated from supplied abundances"))
    }
  }
  if (missing(x) & !missing(size) & !missing(rich)){
    S <- rich
    N <- size
  }
  if (missing(x) & missing(size) & missing(rich)){
    stop("Please provide size and species number or a vector of species abundances")
  }
  if (missing(start.value)){
    f1 <- function(a) {
      S + a*log((a/(a + N)))
    }
    sol <- uniroot(f1, interval = c(1/N, N))
    alfa <- sol$root
    X <- N/(N + alfa)
  } else{
    alfa <- start.value
  }
  if (!missing(x)){
    if (missing(trunc)){
      LL <- function(alpha) -sum(dls(x, N, alpha, log = T))
    } else{
      LL <- function(alpha) -sum(trunc("dls", x, N, alpha, trunc = trunc, log = T))
    }
    result <- mle2(LL, start = list(alpha = alfa), data = list(x = x), method = "Brent", lower = 0, upper = upper, ...)
    if(abs(as.numeric(result@coef) - upper) < 0.0000001) warning("check the upper limit of the function")
    new("fitsad", result, sad = "ls", trunc = ifelse(missing(trunc), NaN, trunc))
  }
  else new("fitsad", coef = c(alpha = alfa), fullcoef = c(alpha = alfa), sad = "ls", trunc = ifelse(missing(trunc), NaN, trunc))
}