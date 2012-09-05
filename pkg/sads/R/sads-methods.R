setGeneric("points")
setMethod("plot", "rad",
          function(x, ...){
            dots <- list(...)
            dots$log <- "y"
            if(!"xlab" %in% names(dots)) dots$xlab = "Species Rank"
            if(!"ylab" %in% names(dots)) dots$ylab = "Species Abundance"
            do.call(plot, c(list(x = x[, 1], y = x[, 2]), dots)) 
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
            x.hist <- rep(as.integer(as.character(x$octave)), as.integer(as.character(x$Freq)))
            hist(x.hist, col = "gray", main = "", ylab = "Frequency", xlab = "Abundance class", breaks = c((min(as.integer(as.character(x$octave)))-1), as.integer(as.character(x$octave))), ...)
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
          function(object, which=c("octaves","rad"), ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            y <- object@data$x
            cf <- as.list(object@coef)
            oct.df <- octav(y)
            rad.df <- rad(y)
            oct.pred <- octavpredt(y, sad = object@sad, coef = cf, ...)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            rad.pred <- radpredt(y, sad=object@sad, coef=cf, ...)
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
              points(rad.pred, ...)
            }
          }
          )
