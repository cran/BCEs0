\name{data}
\alias{data}
\alias{acupuncture}
\docType{data}
\title{
Bayesian model for the cost-effectiveness of acupuncture
}
\description{
This data set contains a sample of the data used in the Bayesian analysis used to 
model the cost-effectiveness of acupuncture in the management of chronic headache.
Data consists in the measure of effectiveness (in terms of annual QALYs) for the
two arms of the trial ($t=0$ indicates usual care, while $t=1$ indicates usual care
and in addition up to 12 acupuncture treatments over three months from appropriately 
trained physiotherapists). Also, the dataset contains the values of the upper bounds
for the average and standard deviation of the cost distribution. 
}
\usage{data(acupuncture)}
\format{
  The format is:
 chr "acupuncture"
}
\source{
Adapted from Wonderling, D., A. Vickers, R. Grieve, and R. McCarney (2004). 
Cost effectiveness analysis of a randomised trial of acupuncture for chronic
headache in primary care. British Medical Journal 328 (7442), 747-752
}
\references{
Baio G. (2013). Bayesian models for cost-effectiveness analysis in the presence of 
structural zero costs. http://arxiv.org/pdf/1307.5243v1.pdf
}
\seealso{
\code{\link{bces0}}
}
\examples{
data(acupuncture)
}
\keyword{datasets}
