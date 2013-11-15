\name{writeModel}
\alias{writeModel}
\title{
writeModel
}
\description{
This function selects suitable bits of JAGS code to build the model file encoding the
selected distributional assumptions for the cost and effectiveness variables (and for 
the selection model) 
}
\usage{
writeModel(dist.c, dist.e, dist.d, model.file)
}
\arguments{
  \item{dist.c}{
A text string defining the selected distribution for the costs. Available options are 
Gamma ("gamma"), log-Normal ("logn") and Normal ("norm")
}
  \item{dist.e}{
A text string defining the selected distribution for the measure of effectiveness. 
Available options are Beta ("beta"), Gamma ("gamma"), Bernoulli ("bern") and Normal 
("norm")
}
  \item{dist.d}{
A text string defining the selection model. Possible choices are "cov.cauchy" or 
"cov.norm" (used when individual covariates are available and can be used to estimate 
the probability of zero costs) and "int" (when no covariate is available and an 
intercept-only model is fitted). The function writes a text file in the current 
working directory, including the relevant bits of code, that can be then passed to 
the call to the function jags to run the MCMC simulations in background 
}
  \item{model.file}{
A string with the name of the model file to which the JAGS code is saved
}
}
\value{
Writes out the file with the selected distributional assumptions to the file 
"model.txt" in the current working directory
}

\references{
Baio G. (2013). Bayesian models for cost-effectiveness analysis in the presence 
of structural zero costs. http://arxiv.org/pdf/1307.5243v1.pdf
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{bces0}}
}
\examples{
## To be added here
}
\keyword{JAGS}
\keyword{Bayesian models}% __ONLY ONE__ keyword per line
