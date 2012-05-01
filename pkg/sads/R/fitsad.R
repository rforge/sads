fitsad <- function(x,sads=c("ls","poilog"),...){
  if(any(sads=="ls")) ls.m <- fitlogser(x,...)
  if(any(sads=="poilog")) p.m <- fitpoilog(x,...)
  new("fitsadlist",list(logseries=ls.m,poilog=p.m))
}

