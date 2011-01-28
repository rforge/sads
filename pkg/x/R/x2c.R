x2c <- function(x) {
  .C("quadrado", as.integer(x), PACKAGE="x")
}

.First.lib <- function(libs, x) {
  library.dynam("x", x, libs)
}