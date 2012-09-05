Power, Zipf and Zipf-Mandelbrot distribution
x <- 1:10
#q = 0  N finite Mandelbrot -> Zipf
dzipf(x, length(x), 2)
dmand(x, length(x), 2, 0)

#q = 0  N infinite Mandelbrot -> Power or Zeta distribution
dpower(x, 2)
dmand(x, 1000000, 2, 0)


dmand(x, length(x), 2, 0)
pmand(x, length(x), 2, 0)
pzipf(x, length(x), 2)

dmand(x, 1000000, 2, 0)
pmand(x, 1000000, 2, 0)
ppower(x, 2)

dmand(x, length(x), 2, 1)
pmand(x, length(x), 2, 1)

dmand(x, 1000000, 2, 1)
pmand(x, 1000000, 2, 1)

p <- pmand(x, length(x), 2, 1)
qmand(p, length(x), 2, 1)

1p <- pzipf(x, length(x), 2)
qzipf(p, length(x), 2)

x<-1.5:10.5
dpareto(x, 0.8091838, 1.5)
VGAM::dpareto(x, 1.5, 0.8091838)
ppareto(x, 0.8091838, 1.5)
VGAM::ppareto(x, 1.5, 0.8091838)
p <- ppareto(x, 0.8091838, 1.5)
VGAM::qpareto(p+0.000001, 1.5, 0.8091838)
p <- ppareto(x, 0.8091838, 1.5, lower = F, log = T)
qpareto(p, 0.8091838, 1.5, lower = F, log = T)

fitmand(x1, trunc = 1)

fitpareto(x1, trunc = 1)
fitnbinom(samp1)
fitnbinom(samp2)

####FITSADS
#Discretas
fitgeom(samp1)
fitsads(samp1, "geom")
fitgeom(samp1, trunc = 0)
fitsads(samp1, "geom", trunc = 0)

fitls(samp1)
fitsads(samp1, "ls")
fitls(x1, trunc = 1)
fitsads(x1, "ls", trunc = 1)

#Continuas
fitgamma(samp1)
fitsads(samp1, "gamma")
fitgamma(samp1, trunc = 0)
fitsads(samp1, "gamma", trunc = 0)
fitgamma(x1, trunc = 1)
fitsads(x1, "gamma", trunc = 1)


####QQSAD
## Poilog
qqsad(samp1.pln)
qqsad2(samp1.pln)
## Logserie
qqsad(samp1.ls)
qqsad2(samp1.ls)
## Gamma
qqsad(samp1.gm)
qqsad2(samp1.gm)
## Power
qqsad(samp1.pw)
qqsad2(samp1.pw)

####PPSAD
## Poilog
ppsad(samp1.pln)
## Logserie
ppsad(samp1.ls)
## Gamma
ppsad(samp1.gm)
## Power
ppsad(samp1.pw)

