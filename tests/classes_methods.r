## Criando classes e metodos para o pacote sads

setClass("fitsad",representation("mle2",sad="character"))
setClass("fitsadlist",representation("list"))
setClass("octav", representation("data.frame"))
setClass("rad", representation("data.frame"))
setMethod("plot","rad",plot.rad)
setMethod("points","rad",points.rad)
setMethod("plot","octav",plot.octav)
setMethod("points","octav",points.octav)
setMethod("plot","fitsad",plot.fitsad)
