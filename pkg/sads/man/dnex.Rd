\name{dnex}
\alias{dnex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compound Negative Binomial-Exponential Distribution
}
\description{
Density for the Negative Binomial-Exponential probability distribution.
}
\usage{
dnex(y, frac, rate, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector of (non-negative integer) quantiles.
}
\item{frac}{
  single numeric '0<frac<=1'; fraction of the population or community sampled (see details)
}
\item{rate}{
  Vector of (non-negative) rates of the exponential distribution of the
  sampled population (see details).
}
\item{k}{
	  dispersion parameter of the negative binomial. Must be strictly positive, need not be
          integer.
}
}

\details{
A compound Negative Binomial-exponential distribution is a Negative Binomial probability distribution
where its expected value is frac*n, where n is a random variable with exponential distribution. The density function
is given by Green & Plotkin (2007)

In ecology, this distribution gives the probability that a species has
an abundance of x individuals in a negative binomial sample of a fraction 'frac' of
the community. In the community the species abundances are independent
random variables that follow an exponential density function. Using of a
negative binomial is justified if individuals of each species are
aggregated across sampling units. The degree of aggregation depends on
the argument 'k', which is the dispersion parameter of the negative
binomial (argument 'size' in dnbinom)

Hence, a Negative Binomial-exponential distribution is a model for species
abundances distributions (SAD) under the assumptions that (a) species
abundances in the community are independent identically distributed
exponential variables, (b) sampling is a Negative Binomial proccess with expected
value frac*n, (c) the sampling is with replacement, or the fraction
sampled is small enough to approximate a sample with replacement.

}
\value{
(log) density of the (zero-truncated) density
}
\references{

  Pielou, E.C. 1977. \emph{Mathematical Ecology}. New York: John Wiley
  and Sons.
  
  Green,J. and Plotkin, J.B. 2007 A statistical theory for sampling species abundances. \emph{Ecology Letters 10}:1037--1045
}
\author{
Paulo I. Prado & Cristiano Strieder
}
\note{
z
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
dexp, dpois, dnbinom
}
