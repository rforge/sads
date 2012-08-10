library(untb)
library(MASS)
library(vegan)
## Result files from TeTame
(tet.ex.out <- read.table("thetame_test_out.txt", header=T))
## Fit of test data of TeTame by untb package ##

## Data: Example datasets of TeTame (zeroes ommited)
tet.ex1 <- count(c(1,8,80,2,1,20))
tet.ex2 <- count(c(15,8,8,20,10,2))
tet.ex3 <- count(c(15,81,80,2,1,2))
## 1. Ewens sampling formula
## Thetas mles : all match those calculated with TeTame
optimal.theta(tet.ex1)
optimal.theta(tet.ex2)
optimal.theta(tet.ex3)
## Log-likelihoods: none matches
theta.likelihood(optimal.theta(tet.ex1), x=tet.ex1)
theta.likelihood(optimal.theta(tet.ex2), x=tet.ex2)
theta.likelihood(optimal.theta(tet.ex3), x=tet.ex3)

## 2. Etienne's sampling formula
##logkda calculations
L1 <- logkda.R(tet.ex1, use.brob=TRUE)
L2 <- logkda.R(tet.ex2, use.brob=TRUE)
L3 <- logkda.R(tet.ex3, use.brob=TRUE)
## Estimation of m and theta: good matching to estimates from TeTame
(tetam1.zsm <- optimal.params(tet.ex1,log.kda=L1,give=FALSE))
(tetam2.zsm <- optimal.params(tet.ex2,log.kda=L2,give=FALSE))
(tetam3.zsm <- optimal.params(tet.ex3,log.kda=L3,give=FALSE)) 
## LogLikelihood: do not match with an order of magnitude
etienne(theta=tetam1.zsm[1],m=tetam1.zsm[2],D=tet.ex1,log.kda=L1)
etienne(theta=tetam2.zsm[1],m=tetam2.zsm[2],D=tet.ex2,log.kda=L2)
etienne(theta=tetam3.zsm[1],m=tetam3.zsm[2],D=tet.ex3,log.kda=L3)
## Comparing with likelihood for a lognormal
logLik(fitdistr(tet.ex1, "lognormal")) # magnitude compatible with log-likelihood values form Tetame
logLik(fitdistr(tet.ex2, "lognormal"))
logLik(fitdistr(tet.ex3, "lognormal"))
