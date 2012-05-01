plot.fitsad <- function(x,which=c("octaves","rad"), ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
  y <- x@data$x
  cf <- as.list(coef(x))
  oct.df <- octav(y)
  rad.df <- rad(y)
  oct.pred <- octavpred(y,sad=x@sad,coef=cf,...)
  oct.ymax <- max(c(oct.df[,3],oct.pred[,3]),na.rm=TRUE)
  rad.pred <- radpred(y,sad=x@sad,coef=cf,...)
  rad.ylim <- range(c(rad.df[,2],rad.pred[,2]),na.rm=TRUE)
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if("octaves"%in%which){
    plot(oct.df,ylim=c(0,oct.ymax),...)
    points(oct.pred,...)
  }
  if("rad"%in%which){
    plot(rad.df,ylim=rad.ylim,...)
    points(rad.pred,...)
  }
}
