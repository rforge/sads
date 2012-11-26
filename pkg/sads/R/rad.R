rad <- function(x){
  if(is(x, "fitsad"))
    new("rad", data.frame(rank=1:length(x@data$x), abund=sort(x@data$x, decreasing=T)))
  if(is(x,"fitrad"))
    x@rad.tab
  else if(is(x,"numeric"))
    new("rad", data.frame(rank=1:length(x), abund=sort(x, decreasing=T)))
}
