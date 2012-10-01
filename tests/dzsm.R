system("R CMD SHLIB zsm.c")
dyn.load("zsm.so")

sn <- function(n, J, m, theta, precisao = 0.01){
  lengn <- length(n)
  res <- rep(0, lengn)
  fun <- .C("sn", as.integer(n), as.integer(lengn), as.integer(J), as.double(m),
            as.double(theta), as.double(precisao), result = res)
  fun[["result"]]
}

system.time((teste <- alonso10(1:150, J=150, m=0.05, theta=2)))
sum(alonso10(1:150, 150, 0.05, 2))

system.time((teste2 <- sn(1:150, 150, 0.05, 2)))
sum(teste2)
plot(teste, teste2)
abline(0,1)

alonso14(150, 150, 0.05, 2)
alonso14A(150, 150, 0.05, 2)

(alonso14(n=1:5, J=15000, m=0.77, theta=41))
(alonso14A(n=1:5, J=15000, m=0.77, theta=41))