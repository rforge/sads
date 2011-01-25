x2c <- function(x) {
  .C("quadrado", as.integer(x), PACKAGE="sads")
}

.First.lib <- function(libs, sads) {
  library.dynam("sads", sads, libs)
}