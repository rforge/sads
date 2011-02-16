\name{depois}
\alias{depois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probability density function for the Poisson-Exponencial combination.
}
\description{
This function describes the probability density as a combination of two parts. 
An exponencial distribution describing the species abundance in the comunity and 
a Poisson distribution describing the sample proccess.
}
\usage{
depois(x, a = 0.1, lambda = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Number of individuals in the community.
}
  \item{a}{
Fraction of the community sampled.
}
  \item{lambda}{
Exponencial parameter. Specie dominance in the comunity.
}
}
\details{
z
}
\value{
z
}
\references{
Prado, P. I. K. L., Species Abundance Distribution, 
Analitical advances in the understanding of an ecological basic pattern.
}
\author{
Cristiano
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
function(x, a = 0.1, lambda = 0.5) {
  
  b <- a^(x-1)
  m <- lambda*b
  n <- (lambda + a)^x
  m/n
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
