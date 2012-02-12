\name{dsad}
\alias{dsad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Numeric approximation to Compound Distributions from Poisson or Negative Binomial Sampling}
\description{Density probability of a Poisson (or Negative Binomial) distribution with
  expected value (parameters 'lambda' or 'mu') distributed as a continuous random variable

}
\usage{
dsad(y, frac, sad, samp="Poisson", k=0.5, log=FALSE, upper=0.9999, trunc=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{
  vector of (non-negative integer) quantiles. For  modelling species abundance
  distributions (SADs), the values are the counts of individuals of each species in the sample
}
\item{frac}{
single numeric '0<frac<=1'; fraction of the community sampled
}

\item{sad}{
name; root name of community sad distribution - e.g., 'lnorm' for the lognormal distribution; 'geom' for the geometric distribution.  
}

\item{samp}{
  a character string indicating the sampling process (Poisson or
  negative binomial). One of '"pois"'
  (default) or '"nbinom"', can be abbreviated. Currently only Poisson
  sampling is implemented.
}

\item{k}{
positive; size parameter for the sampling binomial negative. To be used
when binomial negative sampling is implemented. 
}

\item{log}{
  logical; if TRUE, probabilities p are given as log(p).
}

\item{upper}{
proportion of the community species abundance distribution used to define the upper limit of the integration. '0>upper<=1',  see details.
}

\item{trunc}{
  non-negative integer; truncation point of the compund distribution (see details). 
}

\item{...}{
parameters to be passed to the probability function defined
by the argument 'sad'
}
}
\details{
This function returns the density probability of a Poisson or Negative
Binomial distribution that have their expected value distributed according another
continuous probability density. These expected values correspond to the parameters 'lambda'
for the Poisson and 'mu' for the negative binomial. The resulting
compound distribution is discrete.

In ecology these compound distributions are used to describe
the abundance distributions of species (SADs) in samples taken from biological communities.
The expected species abundances in the sampled community follow a
probability distribution given by the argument 'sad'. A fraction 'frac' of
the community is sampled with replacement, thus the expected abundance in the sample of each
species is frac*n, where n is the species' expected abundance in the
community.

Two sampling processes can be used: a Poisson sample,
where individuals are sampled independently, and Negative Binomial
sample, where individuals are aggregated over sampling units. Currently
only the Poisson sampling is implemented.

In general terms, this function return Poisson or Negative
Binomial probabilities densities where the parameter that correspond to
the mean (= expected value of the variable) follows a probability
distribution given by the argument 'sad'.  Some of these compound
distributions are well-known models for sads (species abundance distributions), like the poisson-lognormal, or the
poisson-gamma distributions. The poisson-gamma compounding results in
the negative binomial distribution, and the Fisher log-series is a limiting case of this distribution.

Taking g(n) as the distribution of abundances of species in the sampled community and
f(y|n) as the conditional distribution of abundances in the sample given an
abundance n in the community, the
general expression for the probability distribution of abundances in the
sample is

f(y) = integral_0^Inf g(n) f(y|n) dn

(Pielou 1977, Green & Plotkin 2007). Function 'dsad' uses numerical integration through function 'integrate' to solve this
expression. To avoid errors for large values of 'y', the upper limit of
integration is substituted by the quantile of g(n) at the probability
level given by the argument 'upper'. The default value is 0.9999, thus the upper limit of
the numerical integration correspond to the abundance value that
encompasses 99.99\% of the species in the community. Although for many
integration problems the use of finite limit values should be avoided
(see help for 'integrate'), in the present case this procedure resulted in
best approximation to analytical values (see examples below)

These compound probability distributions can be used as SAD models under
three key assumptions: (a) species abundances in the community are independent identically distributed
variables that follow a theoretical distribution given by the argument
'sad'; (b) sampling is a Poisson or negative binomial proccess with expected
value frac*n, where n are the abundances in the community (c); individuals are sampled with replacement, or the
fraction of total individuals sampled is small enough to approximate a sample with replacement. See
Engen (1977) and Alonso et al. (2008) for critic evaluations.

Setting a value to 'trunc' truncates the probability distribution at
this value, that is probabilities will be calculated as
f(y)/(1-f(y)). For species-abundance distributions the most used option is to truncate at zero, since species with zero abundance usually are not known.
}
\value{
vector of (log) probabilities
}
\references{
  
  Alonso, D. and Ostling, A., and Etienne, R. S. 2008 The implicit
  assumption of symmetry and the species abundance
  distribution. \emph{Ecology Letters, 11}: 93-105.

  Engen, S. 1977. Comments on two different approaches to the analysis
  of species frequency data. \emph{Biometrics, 33}: 205-213.

  Pielou, E.C. 1977. \emph{Mathematical Ecology}. New York: John Wiley
  and Sons.
  
  Green,J. and Plotkin, J.B. 2007. A statistical theory for sampling
  species abundances. \emph{Ecology Letters 10}:1037--1045

}
\author{
Paulo I. Prado & Cristiano Strieder
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  dpoix, dpoig, ... for densities of compound distributions for SADs
  using analytical given by Green & Plotkin (2007)
  dpois, dnbinom, integrate
}

\examples{

##Poisson-exponential distribution
## Community sad is a exponential with mean abundance = 1000 individuals
## 5\% of the community sampled
## Probability for abundances 0 to 50
x <- 0:50
y1 <- dsad(x,frac=0.05,exp,samp="Poisson",rate=1/1000)
## Comparing with the analytical Poisson-exponential
y2 <- dpoix(x, frac=0.05, rate=1/1000)
plot(x,y1, xlab="y", ylab="prob", ylim=range(c(y1,y2)), cex=0.8)
points(x,y2,pch=3,cex=0.5, type="b")
legend("topright", c("analytical","numeric"), pch=c(1,3))
}
