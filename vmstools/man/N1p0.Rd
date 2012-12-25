\name{N1p0}
\alias{N1p0}
\title{2-dimensional Gaussian distribution with two-peak possibility}
\description{
Computes the 2-dimensional Gaussian distribution depending on the distance from the mean, and
standard deviation with adjustable mean and shift from mean into two peaks instead of 1
}
\usage{
N1p0(x, mu, sig, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Observed value, deviance from mean}
  \item{mu}{Mean value}
  \item{sig}{Standard deviation of value}
  \item{p}{Shift from mean into two peaks}
}
\references{EU lot 2 project}
\author{Niels T. Hintzen}
\note{Function is called inside plotCIinterpolation()}
\examples{

x   <- matrix(c(5,4,3,4,5,4,3,2,3,4,3,2,1,2,3,
                3,2,1,2,3,4,3,2,3,4,5,4,3,4,5),nrow=6,ncol=5,byrow=TRUE)
mu  <- 0
sig <- x^0.5
p   <- 0

res <- matrix(N1p0(x,mu,sig,p),nrow=6,ncol=5)

  #Plot the results
image(res)
}

