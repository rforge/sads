## Funcao names.args V& R pp.46
name.args <- function(...){
  dots <- as.list(substitute(list(...)))[-1]
  nm <- names(dots)
  fixup <- if(is.null(nm)) seq(along=dots) else nm==""
  dep <- sapply(dots[fixup],function(x)deparse(x)[1])
  if(is.null(nm))dep
  else{nm[fixup] <- dep;nm}
}


