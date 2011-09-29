\name{dneg}
\alias{dneg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compound Negative Binomial-Gamma Distribution
}
\description{
Density for the Negative Binomial-Gamma compound probability distribution.
}
\usage{
dneg(y, frac, rate, shape, k)
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
  vector of (non-negative) rates of the gamma distribution of the
  sampled community or population (see details). Inverse of the scale
  parameter of gamma. Must be stricly positive.
}

\item{shape}{
  the shape parameter of the gamma distribution of the
  sampled population (see details).  Must be positive.
}

\item{k}{
	  dispersion parameter of the negative binomial. Must be strictly positive, need not be
          integer.
	}
	
}
\details{
A compound Negative Binomial-gamma distribution is a Negative Binomial probability distribution
where its expected value is frac*n, where n is a random variable with gamma distribution. The density function
is given by Green & Plotkin (2007)

In ecology, this distribution gives the probability that a species has
an abundance of x individuals in a negative binomial sample of a fraction 'frac' of
the community. In the community the species abundances are independent
random variables that follow an gamma density function. Using of a
negative binomial is justified if individuals of each species are
aggregated across sampling units. The degree of aggregation depends on
the argument 'k', which is the dispersion parameter of the negative
binomial (argument 'size' in dnbinom)

Hence, a Negative Binomial-gamma distribution is a model for species
abundances distributions (SAD) under the assumptions that (a) species
abundances in the community are independent identically distributed
gamma variables, (b) sampling is a Negative Binomial proccess with expected
value frac*n, (c) the sampling is with replacement, or the fraction
sampled is small enough to approximate a sample with replacement.
}

\value{
(log) density of the (zero-truncated) density
}
\references{
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
dgamma, dnbinom
}
