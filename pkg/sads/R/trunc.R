trunc<-function(f, x, trunc, log = FALSE, ...){
  #dots <-c(list(...))
  first<-strsplit(f, NULL)[[1]][1]
  fun<- paste(strsplit(f, NULL)[[1]][-1], collapse="")
  switch(first, 
         d = dtrunc(fun, x, trunc, log, ...), 
         p = ptrunc(fun, x, trunc, ...),
         q = qtrunc(fun, x, trunc, ...))
}