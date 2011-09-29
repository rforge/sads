\name{dpoig}
\alias{dpoig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compound Poisson-Exponential Distribution
}
\description{
Density for the Poisson-gamma compound probability distribution.
}
\usage{
dpoig(y, frac, rate, shape, trunc=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Vector of (non-negative integer) quantiles.
}
\item{frac}{
  Probability of success in each trial. '0 < prob <= 1'
}
\item{rate}{
  Vector of (non-negative) rates of the exponential distribution of the
  sampled population (see details).
}
\item{shape}{
	  The shape parameter of the gamma mixing
          distribution.  Must be strictly positive, need not be
          integer.
}
\item{trunc}{
  The most used option is to truncate at zero, since species with zero abundance usually are not known.
}
}
\details{
A compound Poisson-exponential distribution is a Poisson probability distribution
where its single parameter, the process mean rate, is frac*n, where n
is a random variable with exponential distribution. The density function
is

p(x) = rate*frac^x / (frac + rate)^(x+1)

In ecology, this distribution gives the probability that a species has
an abundance of x individuals in a random sample of a fraction frac of
the community. In the community the species abundances are independent
random variables that follow an exponential density function.

Hence, a Poisson-exponential distribution is a model for species
abundances distributions (SAD) under the assumptions that (a) species
abundances in the community are independent identically distributed
exponential variables, (b) sampling is a Poisson proccess with expected
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

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
dexp, dpois
}
