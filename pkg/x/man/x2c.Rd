\name{x2c}
\alias{x2c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Square made by C
}
\description{
z
}
\usage{
x2c(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
z
}
\value{
z
}
\references{
z
}
\author{
Strieder
}
\note{
z
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
z
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x) .C("quadrado", as.integer(x), PACKAGE="x")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
