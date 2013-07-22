\name{bces0}
\alias{bces0}
\alias{bces0.default}
\title{
Bayesian Cost-Effectiveness models in the presence of structural zeros
}
\description{
Writes a model file encoding the distributional assumptions, calls JAGS in background 
and perform the Bayesian analysis of the selected model. 
}
\usage{
bces0(data, dist.c = c("gamma", "logn","norm"), 
      dist.e = c("beta", "gamma", "bern", "norm"), 
      n.iter = 10000, n.burnin = 5000, n.chains = 2) 

\method{bces0}{default}(data, dist.c = c("gamma", "logn","norm"), 
      dist.e = c("beta", "gamma", "bern", "norm"), 
      n.iter = 10000, n.burnin = 5000, n.chains = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A named list including values for the variables e0 (measure of effectiveness for the 
subjects in treatment arm $t=0$), e1 (effectiveness for the subjects in $t=1$), c0 
(individual costs in $t=0$), c1 (individual costs in $t=1$), H.psi (upper bound for the 
prior distribution of the mean of the positive costs), H.zeta (upper bound for the 
prior distribution of the standard deviation of the positive costs). Additional 
optional elements are X0 (a matrix of covariates for $t=0$) and X1 (a matrix of 
covariates for $t=1$) that can be used to estimate the selection model for null costs 
}
  \item{dist.c}{
A text string defining the selected distribution for the costs. Available options are 
Gamma ("gamma"), log-Normal ("logn") and Normal ("norm")
}
  \item{dist.e}{
A text string defining the selected distribution for the measure of effectiveness. 
Available options are Beta ("beta"), Gamma ("gamma"), Bernoulli ("bern") and Normal
("norm")
}
  \item{n.iter}{
Number of iterations to be run in JAGS (default = 10000)
}
  \item{n.burnin}{
Number of iterations to be used as burn-in for the MCMC procedure (default = 5000)
}
  \item{n.chains}{
Number of Markov chains to be run (default = 2)
}
}
\value{
An object containing the following elements
\item{mod}{A "rjags" objects with the results of the MCMC simulations run using JAGS}
\item{params}{A vector including the parameters being monitored}
\item{dataJags}{A list contaning the data needed to run the MCMC simulations}
\item{inits}{A function used to initialise the random nodes in the model}
}
\references{
Baio G. (2013). Bayesian models for cost-effectiveness analysis in the presence of 
structural zero costs. http://arxiv.org/pdf/1307.5243v1.pdf
}
\author{
Gianluca Baio
}
\examples{
data(acupuncture)
m <- bces0(data,dist.c="gamma",dist.e="beta",n.iter=1000,n.burnin=500,n.chains=2)
print(m)
plot(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{JAGS}
\keyword{Markov Chain Monte Carlo}
\keyword{Bayesian models for cost-effectiveness analysis}
