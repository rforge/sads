setMethod("plot","rad",
          function(x,...){
            dots <- list(...)
            dots$log="y"
            if(!"xlab"%in%names(dots)) dots$xlab="Species Rank"
            if(!"ylab"%in%names(dots)) dots$ylab="Species Abundance"
            do.call(plot,c(list(x=x[,1],y=x[,2]),dots)) 
          }
          )

setMethod("points","rad",
          function(x,...){
            dots <- list(...)
            if(!"type"%in%names(dots)) dots$type="l"
            if(!"col"%in%names(dots)) dots$col="blue"
            do.call(points,c(list(x=x[,1],y=x[,2]),dots)) 
          }
          )

setMethod("plot","octav",
          function(x,...){
            barplot(height=x$Freq, space=0, ylab="Frequency", xlab="Abundance class",...)
            axis(1, at=x$octave, labels=x$upper)
          }
          )

setMethod("points","octav",
          function(x,...){
            dots <- list(...)
            if(!"type"%in%names(dots)) dots$type="b"
            if(!"col"%in%names(dots)) dots$col="blue"
            X <- c(0,as.integer(as.character(x$octave)))
            X <- X[-length(X)]+diff(X)/2
            do.call(points,c(list(x=X,y=x$Freq),dots))
          }
          )

setMethod("plot","fitsad",
          function(x,which=c("octaves","rad"), ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
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
          )
