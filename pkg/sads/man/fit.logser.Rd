\name{fit.logser}
\alias{fit.logser}
\title{Fits Fisher's Logseries Distribution}
\description{
Given a vector of species abundances or the total number of
species and individuals in a sample, fits Fisher's logseries distribution
by maximum likelihood.
}
\usage{
fit.logser( x, size , rich )
}

\arguments{
  \item{x}{
    Vector of (non-negative integer) abundances of species in a sample.
  }
  \item{size}{
    Total number of individuals in the sample. Ignored if x is supplied.
  }
  \item{rich}{
    Total number of species in the sample. Ignored if x is supplied
  }
}

\details{

  The Fisher logseries is a limiting case of the Negative Binomial where
  the dispersion parameter of the negative binomial tends to zero. It
  was originally proposed by Fisher (1943) to relate the expected number
  of species in a sample from a biological community to the sample size as:

  S = alpha * log(1 + N/alpha)

  Where alpha is the single parameter of the logseries distribution,
  often used as a diversity index. From this relation follows that the
  expected number of species with x individuals in the sample is
  
  S(x) = alpha*X^x/x

  Where X is a function of alpha and N, that tends to one as the sample
  size N increases:

  X = N / (alpha+ N)
  
  This function implements the original numerical recipe by Fisher to
  fit the logseries distribution, given a vector of abundances, or ,
  alternatively, the number of individuals and species in the sample. Alonso et al. (2008,
  supplementary material) showed that this recipe gives the maximum
  likelihood estimate of Fisher's alpha. The function also
  calculates the log-likelihood and Akaike information index for the
  fit.

  Fitting is done through numerical optmization with the 'uniroot'
  function, following the code of the function 'fishers.alpha' of the
  untb package. After that, the fitted value of the alpha parameter is
  use to get the Log-likelihood the logseries density
  function 'dls'. 
}

\value{
  A named list with components
  \item{summary}{
    named list with components 'N' and 'S' that give the number of
    individuals and species in the sample; 'alpha' and 'X' give the estimated
    value of Fisher's alpha and X parameter (see details); 'LogLik' and
    'AIC' give the log-likelihood and Akaike Information Index for the
    fitted logseries.
  }
  \item{sol}{
    list returned by the function 'uniroot'
  }
}

\references{

  Alonso, D. and Ostling, A., and Etienne, R. S. 2008 The implicit
  assumption of symmetry and the species abundance
  distribution. \emph{Ecology Letters, 11}: 93-105.
  
  Pielou, E.C. 1977. \emph{Mathematical Ecology}. New York: John Wiley
  and Sons.

  Fisher, R.A, Corbert, A.S. and Williams, C.B. (1943) The Relation
  between the number of species and the number of individuals in a
  random sample of an animal population. \emph{The Journal of Animal Ecology, 12(1)}: 42--58.

  Green,J. and Plotkin, J.B. 2007 A statistical theory for sampling
  species abundances. \emph{Ecology Letters 10}:1037--1045
}

\author{
Paulo I. Prado
}


\seealso{
'dpoig', 'dls', 'uniroot'; fishers.alpha in package untb and 'fisherfit'
in package vegan for other functions to fit Fisher's logseries. 
}
