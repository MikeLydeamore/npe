% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npe-calculations.R
\name{calculateNPE}
\alias{calculateNPE}
\title{Calculate value of n-phase exponential curve}
\usage{
calculateNPE(lambda, k, t)
}
\arguments{
\item{lambda}{Vector specifying the initial size of each group}

\item{k}{Vector specifying the decay rates of each group}

\item{t}{Time(s) at which to evaluate curve}
}
\value{
A vector containing the value of the function,
\eqn{\sum_i \lambda_i exp(k_i \times t)}, at each time \eqn{t}.
}
\description{
Calculates the value of an n-phase exponential curve with initial sizes, \eqn{\lambda} and
rate parameters, \eqn{k} at time \eqn{t}. \cr \cr
Note importantly that the rate parameters, \eqn{k}, are \strong{not} negated.
}
\examples{
lambda <- c(300, 200)
k <- c(-0.1, -0.5)
t <- seq(1:10, by=0.1)
calculateNPE(lambda = lambda, k = k, t = t)
}
