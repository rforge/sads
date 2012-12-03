#setGeneric("points")

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
          function(x, which=1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            oct.df <- octav(x)
            rad.df <- rad(x)
            oct.pred <- octavpred(x)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            rad.pred <- radpred(x)
            rad.ylim <- range(c(rad.df[, 2], rad.pred[, 2]), na.rm = TRUE)
            if (ask) {
              oask <- devAskNewPage(TRUE)
              on.exit(devAskNewPage(oask))
            }
            if(1 %in% which){
              plot(oct.df, ylim = c(0, oct.ymax), ...)
              points(oct.pred, ...)
            }
            if(2 %in% which){
              plot(rad.df, ylim = rad.ylim, ...)
              lines(rad.pred, ...)
            }
            if(3 %in% which){
              qqsad(x, ...)
            }
            if(4 %in% which){
              ppsad(x, ...)
            }
          }
          )

setMethod("plot","fitrad",
          function(x, which=1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
            oct.df <- octav(x)
            rad.df <- rad(x)
            oct.pred <- octavpred(x)
            oct.ymax <- max(c(oct.df[, 3], oct.pred[, 3]), na.rm = TRUE)
            rad.pred <- radpred(x)
            rad.ylim <- range(c(rad.df[, 2], rad.pred[, 2]), na.rm = TRUE)
            if (ask) {
              oask <- devAskNewPage(TRUE)
              on.exit(devAskNewPage(oask))
            }
            if(1 %in% which){
              plot(oct.df, ylim = c(0, oct.ymax), ...)
              points(oct.pred, ...)
            }
            if(2 %in% which){
              plot(rad.df, ylim = rad.ylim, ...)
              lines(rad.pred, ...)
            }
            if(3 %in% which){
              qqrad(x, ...)
            }
            if(4 %in% which){
              pprad(x, ...)
            }
          }
          )

## copy of the method in bbmle, with a line added to assure df is not NULL
setMethod("AICc","fitsad",
          function (object, ..., nobs, k = 2){
            L <- list(...)
            if (length(L)) {
              L <- c(list(object), L)
              if (missing(nobs) && is.null(attr(object, "nobs"))) 
                stop("must specify number of observations")
              nobs <- sapply(L, attr, "nobs")
              if (length(unique(nobs)) > 1) 
                stop("nobs different: must have identical data for all objects")
              logLiks <- sapply(L, logLik)
              df <- sapply(L, attr, "df")
              if(is.null(df)) df <- sapply(L, function(object) attr(logLik(object),"df")) ## added to assure that df is not NULL
              val <- -2 * logLiks + k * df * (df + 1)/(nobs - df - 
                                                       1)
              data.frame(AICc = val, df = df)
            }
            else {
              df <- attr(object, "df") 
              if(is.null(df)) df <- attr(logLik(object),"df") ## added to assure that df is not NULL
              c(-2 * logLik(object) + k * df + k * df * (df + 1)/(nobs - 
                                                                  df - 1))
            }
          }
          )
