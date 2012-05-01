radpoints <- function(x,...){
  points(1:length(x),x[order(x,decreasing=T)],...)
}
