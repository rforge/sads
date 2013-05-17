\name{fitsad}
\alias{fitgamma}
\alias{fitgeom}
\alias{fitlnorm}
\alias{fitls}
\alias{fitmzsm}
\alias{fitnbinom}
\alias{fitpareto}
\alias{fitpoilog}
\alias{fitpower}
\alias{fitweibull}
\alias{fitvolkov}
\title{ML fiting of species abundance distributions}
\description{Fits probability distributions for
  abundances of species in a sample or assemblage by maximun likelihood.
}
\usage{
fitsad(x,
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
  \item{sad}{
    character; root name of community sad distribution to be fitted.
    "gamma" for gamma distribution,
    "geom" for geometric distributions (not geometric series rad model,
    \link{\code{dgs}}),
    "lnorm" for logmormal,
    "ls" for Fisher's logseries,
    "mzsm" for Alonso & McKane's neutral metacommunity distribution,
    "nbinom" for negative binomial,
    "pareto" for Pareto distribution,
    "poilog" for Poisson-lognormal distribution,
    "power" for power-law distribution,
    "weibull" for weibull distribution,
    "volkov" for Volkov's et al. neutral community distribution.
 }

    \item{trunc}{
      non-negative integer, trunc > min(x); truncation point to fit a
      truncated distribution.
    }

    \item{start.value}{
      real positive; starting value of free parameters to be passed to \code{mle2}.
  }

  \item{trueLL}{
   logical; use 'true likelihood' for continuous models? (see \link{\code{trueLL}}.
  }

  \item{dec.places}{
    positive integer; number of decimal places used in the measurement of
    the observed values. Observed values will be rounded to this number of
    decimals. This argument defines the measurement precision for
    abundance variables. Used for continuous distributions of
    trueLL=TRUE 
    (see \link{\code{trueLL}}.
  }

  \item{upper}
  real positive; upper bound for the Brent's one-parameter optimization
    method (default), for fits that uses these method by default see
  details and \link{\code{optim}}.
}

\value{
  An object of class 'fitsad' which inherits from 'mle2'
  class and thus has methods for extraction of AIC,
  coeficients and so on (see \code{\link{fitsad-class}})
}

\details{
  \code{fitsad} is only a wraper that calls the functions to fit
  the distribution specified by the argument \code{sad}. Users
  can interchangeably \code{fitsad} or the individuals functions
  detailed below 
  (e.g. \code{fitsad(x, sad="geom", ...)} or \code{fitgeom(x, ...)}.
    
  The distributions are fitted by the
  maximum likelihood method using numerical optimization,
  with \code{mle2}.
  Additionally to methods
  inherited from \code{mle2} objects, the resulting \code{fitsad}
  objects have specific
  methods to plot predicted and observed values (see
  \code{\link{fitsad-class}}.

  Functions \code{fitgamma}, \code{fitlnorm}, \code{fitweibull}, fit the
  standard continuous distributions most used as SADs models.
  Functions \code{power} and \code{pareto} fit power-law continuous distributions
  with one and two-parameters, that have been suggested as SADs models.
  By default, all continuous distributions
  are fit using the
  correct likelihood for independent observations of a
  continuous variable (Lindsey 1999, see \link{\code{trueLL}}).
  Hence, the correct number of decimal places of the abundance
  measure must be given in the argument \code{dec.places}.
  A warning message reminds users about that when \code{trueLL=TRUE}.
  When the abundance measure is the number of individuals of each
  species in the sample,
  the default value \code{dec.places=0}
  is the correct choice. In this case, the probability of a given abundance
  value x will be calculated as CDF(x + 0.5) - CDF( x - 0.5 ), where CDF
  is the cumulative density function of the continuous distribution.

  Functions \code{fitgeom}, \code{fitnbinom}

  Functions \code{fitls} implements the original numerical recipe by Fisher (1943) to
  fit the logseries distribution, given a vector of species abundances.
  Alonso et al. (2008,
  supplementary material) showed that this recipe gives the maximum
  likelihood estimate of Fisher's alpha, the single parameter of the logseries.
  Fitting is done through numerical optimization with the \code{uniroot}
  function, following the code of the function \code{fishers.alpha} of the
  'untb' package. After that, fitted value of the alpha parameter is
  use to get the Log-likelihood from the logseries density
  function \code{dls}, using the function \code{mle2}.

  Function\code{fitbs}
  
  \code{fitmzsm} fits the metacommunity Zero-sum multinomial distribution
  \link{\code{dmzsm}} from the Neutral Theory of Biodiversity.
  Because the total number of individuals (J) is known from the sample,
  the mZSM distribution has a single
  free parameter, theta. By default, \code{fitmzsm} fits mZSM to a vector of abundances
  with Brent's one-dimensional method of optimization (see
  \link{\code{optim}}. The logseries distribution (Fisher et al. 1943)
  is a limiting case of mZSM (Hubbel 2001), and in practice
  the two models provide very similar fits to SADs (see example).

  Function \code{fitvolkov}

}

\references{

  Alonso, D. and McKane, A. J.  2004.  Sampling Hubbell's neutral
  model of biodiversity. \emph{Ecology Letters 7}:901--910

  Alonso, D. and Ostling, A., and Etienne, R. S. 2008 The implicit
  assumption of symmetry and the species abundance
  distribution. Ecology Letters, 11: 93-105.

  Fisher, R.A, Corbert, A.S. and Williams, C.B. (1943) The Relation
  between the number of species and the number of individuals in a
  random sample of an animal population.
  \emph{The Journal of Animal Ecology, 12(1)}: 42--58.
  
  Hubbell, S. P.  2001.  “The Unified Neutral Theory of Biodiversity”.
  Princeton University Press
  
  Lindsey, J.K. 1999. Some stattistical heresies. \emph{The Statistician
    48}(1): 1--40.
}

\source{
  
  }
\seealso{
}

\examples{
}
\author{}
\keyword{}
