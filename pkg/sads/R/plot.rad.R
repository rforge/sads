plot.rad <- function(x,...){
  dots <- list(...)
  dots$log="y"
  if(!"xlab"%in%names(dots)) dots$xlab="Species Rank"
  if(!"ylab"%in%names(dots)) dots$ylab="Species Abundance"
  do.call(plot,c(list(x=x[,1],y=x[,2]),dots)) 
}
