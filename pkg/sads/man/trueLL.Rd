\name{trueLL}
\alias{trueLL}
\title{True likelihood for continuous variables
}
\description{
Calculates the correct likelihood for independent observations of a
continuous variable that follows a given density function.
}
\usage{
trueLL(x, dens, precision=1, log=TRUE,...)
}
\arguments{
  \item{x}{
vector of observed values
}
\item{dens}{
name; root name of the continuous density distribution of the variable -
e.g., 'lnorm' for the lognormal distribution; 'gamma' for the gamma
distribution.
}

\item{precision}{
  positive number; measurement precision of the observed values.
}
\item{log}{
logical; convert probabilities to logs and sum them to get
  log-likelihoods? The log-likelihood value is returned if TRUE and the likelihood
  value if FALSE (quickly tends to zero as number of observation increases).
}
\item{...}{
named arguments to be passed to the probability function defined by the argument 'dens'.
}
}
\details{
The (log)likelihood function is often defined as any function
  proportional to the (sum) product of probability density of the observations. In its original
  formulation, however, the likelihood is proportional to the product of
  probabilities, not probabilities densities (Lindsey 1999, Fisher 1922). For
  continuous variables these probabilities are calculated through
  integration of the probability density from x-D to x+D, where x is the
  observed value and D is half the measurement precision. The likelihood
  function is thus any function proportional to
  
  prod ( integral_(x-D)^(x+D) f(x) dx )

  where D is precision/2, x the value of each observation and f(x) the
  density function. Because products of probabilities quickly tends to
  zero, probabilities are usually converted to their logs and them
  summed to give the log-likelihood function.
}
\value{
(log)-likelihood for the observations, under the probability density
model defined by the argument 'dens' and their parameters, passed
through the argument '...' .
}
\references{
Fisher, R.A. 1922. On the mathematical foundations of theoretical
statistics. \emph{Phil. Trans. R. Soc. Lond. A 222}: 309--368.

Lindsey, J.K. 1999. Some stattistical heresies. \emph{The Statistician
  48}(1): 1--40.
}
\author{
Paulo I. Prado
}


\seealso{
logLik in the package MASS
}

\examples{
##Random sample of a lognormal distribution with precision=0.1
x <- round(rlnorm(500,meanlog=3,sdlog=0.5),1)
## Log-likelihood given by fitdistr
library(MASS)
logLik(fitdistr(x,"lognormal"))
## Which is the sum of log of densities
sum(dlnorm(x,mean(log(x)),sd(log(x)),log=T))
## Correct log-likelihood
trueLL(x,lnorm,precision=0.1,meanlog=mean(log(x)),sdlog=sd(log(x)))
