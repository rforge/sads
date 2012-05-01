rad <- function(x){
  new("rad",
      data.frame(rank=1:length(x),
                 abund=sort(x,decreasing=T))
      )
}
