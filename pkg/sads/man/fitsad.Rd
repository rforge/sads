\name{fitsad}
\alias{fitgamma}
\alias{fitgeom}
\alias{fitlnorm}
\alias{fitls}
\alias{fitnbinom}
\alias(fitpareto}
\alias{fitpoilog}
\alias{fitpower}
\alias{fitweibull}
\alias{fitvolkov}
\title{ML fiting of species abundance distributions}
\description{Fits probability distributions for
  abundances of species in a sample or assemblage
}
\usage{
fitsad(x, sad=
sad =c("gamma","geom","lnorm","ls","mzsm","nbinom","pareto","poilog","power","weibull","volkov"),
trunc, start.value, trueLL = TRUE , dec.places = 0, ...)
fitgamma(x, trunc, start.value, trueLL = TRUE, dec.places = 0, ...  )
fitgeom(x, trunc, start.value, ...  )
fitlnorm(x, trunc, start.value, trueLL = TRUE, dec.places = 0, ...  )
fitls(x, trunc, start.value, upper = length(x), ...  )
fitmzsm(x, trunc, start.value, upper = length(x), ... )
fitnbinom(x, trunc, start.value, ...  )
fitpareto(x, trunc, start.value, trueLL = TRUE, dec.places = 0, upper = 20, ... )
fitpoilog(x, trunc, start.value, ...  )
fitpower(x, trunc, start.value, ...  )
fitweibull(x, trunc, start.value, trueLL = TRUE, dec.places = 0, ...  )
fitvolkov(x, trunc, start.values, ... )
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
\source{
  If starting value is not supplied in \code{fitls}, it is calculated through numerical optimization with the \code{uniroot}
  function, following the code of the function \code{fishers.alpha} of the
  \pkg{untb} package. Since \code{fishers.alpha} implements the
  numerical recipe of Fisher (1943), this starting value is already the
  maximum likelihood estimate for Fisher's alpha Alonso et al. (2008, supplementary material).
  }
\seealso{
}
\examples{
}
\author{}
\keyword{}
