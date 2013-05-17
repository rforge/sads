
\name{fitvolkov}

\alias{fitvolkov}

\title{Fits Neutral Biodiversity Theory distribution}

\description{Fits species abundances distribution
  in a neutral community with
  immigration as deduced
  by Volkov et al 2003 by maximum likelihood.
}

\usage{
fitvolkov( x, trunc, start.values, ... )
}

\arguments{

    \item{x}{
    vector of (non-negative integer) abundances of species in a sample.
  }
  \item{trunc}{
    non-negative integer, trunc > min(x); truncation point to fit a
  truncated distribution.
}
\item{start.value}{
  real positive; starting value for the free parameters theta and m,
  (see \link{\code{dvolkov}}).  
  }
  \item{...}{
    further arguments to be passed to mle2
    }
}

\value{
  \code{fitmzsm} returns an object fo class 'fitsad'.
}

\details{

  There are many proposed models of species abundance distributions (SADs)
  deduced from the Neutral Theory of Biodiversity (Hubbell 2001).
  One of them is the stationary distribution by Volokov et al. (2003),
  for community with neutral drift with immigration. The model has two
  free parameters, the 'fundamental biodiversity number' theta, and the
  immigration rate m (see \link{\code{dvolkov}}).

  \code{fitvolkov} builds on \code{untb::volkov} function to
  fit Volkov's et al SAD model to an vector of abundances,
  with numerical maximum likelihood methods with \code{bbmle::mle2}.
  The fit can be extremely slow even for vectors of moderate size.
  
}

\references{
}

\seealso{
}

\author{Paulo I. Prado and Murilo Dantas Miranda}

\examples{
## birds in an Australian forest
mag5 <- c(103,115,13,2,67,36,51,8,6,61,10,21,
          7,65,4,49,92,37,16,6,23,9,2,6,5,4,
          1,3,1,9,2)
mg5.v <- fitvolkov(mag5)
coef(mg5.v)
par(mfrow=c(2,2))
plot(mg5.v)
par(mfrow=c(1,1))
}
