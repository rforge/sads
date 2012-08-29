trunc<-function(f, x, trunc, coef, ...){
  #dots <-c(list(...))
  first<-strsplit(f, NULL)[[1]][1]
  fun<- paste(strsplit(f, NULL)[[1]][-1], collapse="")
  switch(first, 
         d = dtrunc(fun, x, trunc, coef, ...), 
         p = ptrunc(fun, x, trunc, coef, ...),
         q = qtrunc(fun, x, trunc, coef, ...))
}
