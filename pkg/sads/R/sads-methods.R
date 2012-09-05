setGeneric("points")

setMethod("plot", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"log" %in% names(dots)) dots$log <- "y"
            if(!"xlab" %in% names(dots)) dots$xlab = "Species Rank"
            if(!"ylab" %in% names(dots)) dots$ylab = "Species Abundance"
            if(!"frame.plot" %in% names(dots)) dots$frame.plot = TRUE
            if(!"axes" %in% names(dots)){ 
              do.call(plot, c(list(x = x[, 1], y = x[, 2], axes=FALSE), dots))
              axis(2)
              sc <- axisTicks(range(x[, 1]),nint=10,log=FALSE)
              sc[sc==0] <- 1
              axis(1,at=sc)
            }
            if("axes" %in% names(dots)){ 
              do.call(plot, c(list(x = x[, 1], y = x[, 2]), dots))
            }
            
          }
            )

setMethod("points", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type = "p"
            if(!"col" %in% names(dots)) dots$col = "blue"
            do.call(points, c(list(x = x[, 1], y = x[, 2]), dots)) 
          }
          )

setMethod("lines", "rad",
          function(x, ...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type = "l"
            if(!"col" %in% names(dots)) dots$col = "blue"
            do.call(lines, c(list(x = x[, 1], y = x[, 2]), dots)) 
          }
)

setMethod("plot","octav",
          function(x,...){
            dots <- list(...)
            x.hist <- rep(as.integer(as.character(x$octave)), as.integer(as.character(x$Freq)))
            if(!"col" %in% names(dots)) dots$col = "gray"
            if(!"main" %in% names(dots)) dots$main = ""
            if(!"ylab" %in% names(dots)) dots$ylab = "N of species"
            if(!"xlab" %in% names(dots)) dots$xlab = "Abundance class"
            if(!"axes" %in% names(dots)){ 
              do.call(hist, c(list(x=x.hist,
                   breaks = c((min(as.integer(as.character(x$octave)))-1),as.integer(as.character(x$octave))),
                                   axes=FALSE),dots))
              axis(2)
              n <- as.numeric(as.character(x[,1]))
              axis(1,at=n[seq(1,length(x[,1]),2)],
                   labels=x[seq(1,length(x[,1]),2),2])
            }
            else
              do.call(hist, c(list(x=x.hist,
                   breaks = c((min(as.integer(as.character(x$octave)))-1),as.integer(as.character(x$octave)))),dots))
          }
          )

setMethod("points","octav",
          function(x,...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type="b"
            if(!"col" %in% names(dots)) dots$col="blue"
            X <- c((min(as.integer(as.character(x$octave)))-1), as.integer(as.character(x$octave)))
            X <- X[-length(X)]+diff(X)/2
            do.call(points, c(list(x = X, y = x$Freq), dots))
          }
          )

setMethod("lines","octav",
          function(x,...){
            dots <- list(...)
            if(!"type" %in% names(dots)) dots$type="b"
            if(!"col" %in% names(dots)) dots$col="blue"
            X <- c((min(as.integer(as.character(x$octave)))-1), as.integer(as.character(x$octave)))
            X <- X[-length(X)]+diff(X)/2
            do.call(lines, c(list(x = X, y = x$Freq), dots))
          }
)

setMethod("plot","fitsad",
          function(x, which=c("octaves","rad"), ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            object <- x
            y <- object@data$x
            cf <- as.list(object@coef)
            tr <- object@trunc
            oct.df <- octav(y)
            rad.df <- rad(y)
            if(!is.na(object@trunc))
              oct.pred <- octavpredt(x = y, sad = object@sad, coef = cf, trunc = tr, ...)
            else
              oct.pred <- octavpredt(x = y, sad = object@sad, coef = cf, ...)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            if(!is.na(object@trunc))
              rad.pred <- radpredt(x = y, sad=object@sad, coef=cf, trunc = tr, ...)
            else
              rad.pred <- radpredt(x = y, sad=object@sad, coef=cf, ...)
            rad.ylim <- range(c(rad.df[, 2], rad.pred[, 2]), na.rm = TRUE)
            if (ask) {
              oask <- devAskNewPage(TRUE)
              on.exit(devAskNewPage(oask))
            }
            if("octaves" %in% which){
              plot(oct.df, ylim = c(0, oct.ymax), ...)
              points(oct.pred, ...)
            }
            if("rad" %in% which){
              plot(rad.df, ylim = rad.ylim, ...)
              lines(rad.pred, ...)
            }
            qqsad(object)
            ppsad(object)
          }
          )
