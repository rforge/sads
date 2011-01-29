x2c <- function(x) {
  .C("quadrado", as.integer(x), PACKAGE="tests")
}

.First.lib <- function(libs, tests) {
  library.dynam("tests", tests, libs)
}