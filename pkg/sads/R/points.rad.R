points.rad <- function(x,...){
  dots <- list(...)
  if(!"type"%in%names(dots)) dots$type="l"
  if(!"col"%in%names(dots)) dots$col="blue"
  do.call(points,c(list(x=x[,1],y=x[,2]),dots)) 
}
