\name{fitsad}
\alias{fitgamma}
\alias{fitgeom}
\alias{fitlnorm}
\alias{fitls}
\alias{fitnbinom}
\alias{fitpoilog}
\alias{fitpower}
\alias{fitweibull}
\alias{fitzipf}
\title{ML fiting of species abundance distributions}
\description{Fits basic probability distributions for
  abundances of species in a sample or assemblage
}
\usage{
fitsad(x, distr=
c("gamma","geom","lnorm","ls","nbinom","poilog","power","weibull","zipf"),
trunc, start.value, trueLL, dec.places, ...)
fitgamma(x, trunc, start.value, trueLL, dec.places, ...  )
fitgeom(x, trunc, start.value, ...  )
fitlnorm(x, trunc, start.value, trueLL, dec.places, ...  )
fitls(x, trunc, start.value, ...  )
fitnbinom(x, trunc, start.value, ...  )
fitpoilog(x, trunc, start.value, ...  )
fitpower(x, trunc, start.value, ...  )
fitweibull(x, trunc, start.value, trueLL, dec.places, ...  )
fitzipf(x, trunc, start.value, ...  )
}
\arguments{
  \item{x}{vector of (positive integer) quantiles. In the context of
    sads, some abundance measurement (e.g., number of individuals,
    biomass) of species in a sample or ecological assemblage.
  }
  \item{trunc}{
  }
  \item{start.value}{
  }
  \item{trueLL}{
  }
  \item{dec.places}{
  }
}

\value{
}
\details{
}
\references{
}
\seealso{
}
\examples{
}
\author{}
\keyword{}
