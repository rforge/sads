\name{dpower}
\alias{dpower}
\alias{ppower}
\alias{qpower}
\title{Power distribution
}
\description{
Density, distribution function and quantile function for
power distribution with parameter 's'.
}
\usage{
dpower( x, s, log=FALSE)
ppower( q, s, lower.tail=TRUE, log.p=FALSE)
qpower( p, s, lower.tail= TRUE, log.p=FALSE)
}

\arguments{
  \item{x}{
    vector of (non-negative integer) quantiles. In the context of
    species abundance distributions, this is a vector of abundances of species in a sample.
  }
  \item{q}{vector of (non-negative integer) quantiles. In the context of
    species abundance distributions, a vector of
    abundances of species in a sample.
  }
  \item{p}{
    vector of probabilities.
  }
  \item{s}{
    positive real s > 1; exponent of the power distribution
  }
  \item{log, log.p}{
    logical; if TRUE, probabilities p are given as log(p).
  }
  \item{lower.tail}{
    logical; if TRUE (default), probabilities are P[X <= x],
    otherwise, P[X > x].
  }
}

\details{
  The power density is

  p(x) = x^(-s) / zeta (s)

  Hence p(x) is proportional to a
  negative power of 'x', given by the 's' exponent. The Riemann's zeta
  function is the integration constant.

  The power distribution can be used as species abundance distribution (sad) model, which
  describes the probability of the abundance 'x' of given species in a
  sample or assemblage of species.
 }

\value{
'dpower' gives the (log) density of the density, 'ppower' gives the (log)
distribution function, 'qpower' gives the quantile function.
}

\references{
  Johnson N. L., Kemp, A. W. and Kotz S. (2005) Univariate Discrete
  Distributions, 3rd edition, Hoboken, New Jersey: Wiley. Section 11.2.20. 
}

\author{
Paulo I. Prado and Murilo Dantas Miranda
}

\seealso{
dzeta in VGAM package; rsad for random
generation fitpower for maximum likelihood estimation in the context of species abundance distributions.
}
