points.octav <- function(x,...){
  dots <- list(...)
  if(!"type"%in%names(dots)) dots$type="b"
  if(!"col"%in%names(dots)) dots$col="blue"
  X <- c(0,as.integer(as.character(x$octave)))
  X <- X[-length(X)]+diff(X)/2
  do.call(points,c(list(x=X,y=x$Freq),dots))
}
