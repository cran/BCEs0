\name{plot.bces0}
\alias{plot.bces0}
\title{
Plot method for objects in the class "bces0"
}
\description{
Produces a traceplot for the main nodes in the model
}
\usage{
\method{plot}{bces0}(x, ...)
}
\arguments{
  \item{x}{
The object in the class "BCEs0" obtained by calling the function bces0
}
\item{...}{
Additional arguments affecting the summary produced
}
}
\details{
Returns a traceplot to assess convergence
}
\references{
Baio G. (2013). Bayesian models for cost-effectiveness analysis in the presence of 
structural zero costs. http://arxiv.org/pdf/1307.5243v1.pdf
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{JAGS}
\keyword{Markov Chain Monte Carlo}
\keyword{Bayesian models for cost-effectiveness analysis}
