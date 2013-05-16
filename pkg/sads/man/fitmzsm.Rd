\name{fitmzsm}
\alias{fimzsm}
\title{Fits the metacommunity zero-sum multinomial distribution}
\description{
  Fits Alonso & McKane's mZSM distribution to a vector
  of species abundances by maximum likelihood
}

\usage{
fitmzsm(x, trunc, start.value, upper = length(x), ...)
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
  real positive; starting value for the free parameter theta of the mZSM
  distribution (see \link{\code{dmzsm}}).  
  }
  \item{upper}{
    real positive; upper bound for the Brent's numerical optimization
    method (default), see \link{\code{optim}}.
  }
  \item{...}{
    further arguments to be passed to mle2
    }
}

\details{
  The metacommunity Zero-sum multinomial distribution
  (mZSM)
  describes the probabilities
  of abundances (population sizes)
  in a random sample of size J taken from a
  collection 
  of populations (the metacommunity).
  The total number of individuals in the metacommunity
  is fixed (zero-sum assumption).
  The populations in the metacommunity undergo a stochastic
  birth-death-speciation (or migration) process,
  with the equal demographic rates  (neutrality or ecological
  equivalence assumption, Hubbell 2001).
  \code{fitmzsm} fits the approximation by Alonso et McKane (2004) for mZSM
  for large Poisson samples (J>100). A warning is issued if J<100. 
  
  The mZSM 
  is a model of species
  abundances distribution (SAD)
  in a sample taken from a community under
  the assumptions that (a) species
  abundances in the community follows the stationary distribution
  of a neutral, zero-sum stochastic process of birth, death and speciation (or
  migration);
  (b) sampling is a Poisson proccess with expected
  value aproximated by N(x), (c) individuals are sampled with replacement, or the
  fraction of total individuals sampled is small enough to approximate a
  sample with replacement.

  Because J is known from the sample, the mZSM distribution has a single
  free parameter, theta. By default, \code{fitmzsm} fits mZSM to a vector of abundances
  with Brent's one-dimensional method of optimization (see
  \link{\code{optim}}. 

  The logseries distribution (Fisher et al. 1943)
  is a limiting case of mZSM (Hubbel 2001), and in practice
  the two models provide very similar fits to SADs (see example).
}

\value{
  \code{fitmzsm} returns an object fo class 'fitsad'.
}

\references{
  Alonso, D. and McKane, A. J.  2004.  Sampling Hubbell's neutral
  model of biodiversity. \emph{Ecology Letters 7}:901--910

  Fisher, R.A, Corbert, A.S. and Williams, C.B. (1943) The Relation
  between the number of species and the number of individuals in a
  random sample of an animal population.
  \emph{The Journal of Animal Ecology, 12(1)}: 42--58.
  
  Hubbell, S. P.  2001.  “The Unified Neutral Theory of Biodiversity”.
  Princeton University Press
}

\author{
Paulo I. Prado and Murilo Dantas Miranda
}


\seealso{
  \code{fitsad-class}, \code{dmzsm}, \code{fitls};
  \code{alonso.eqn12}
  in package 'untb' which is based on the exact formulation of mZSM,
  \code{optimal.theta} in package 'untb' for another function for
  estimating theta.
}

\examples{
data(moths) #Fisher's moths data
moths.mzsm <- fitmzsm(moths) # same as fitsad(moths, sad="mzsm")
coef(moths.mzsm) ## Compare with theta=38.9, Alonso&McKanne (2004)

## Diagnostic plots
par(mfrow=c(2,2))
plot(moths.mzsm)
par(mfrow=c(1,1))

## Comparison to logseries
## fit to logseries
moths.ls <- fitsad(moths, sad="ls")
## Model selection
AICtab(moths.ls, moths.mzsm, weights=T)

}
