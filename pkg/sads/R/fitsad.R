fitsad <- function(x, sads=c("ls","poilog"),...){
  if(any(sads == "ls")) ls.m <- fitls(x,...)
  if(any(sads == "poilog")) p.m <- fitpoilog(x,...)
  new("fitsadlist", list(ls = ls.m, poilog = p.m))
}