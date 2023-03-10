% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cumWeightedSum}
\alias{cumWeightedSum}
\title{cumWeightedSum}
\usage{
cumWeightedSum(Index, y, multiplicity, data, K, self_as_NN = FALSE)
}
\arguments{
\item{Index}{The matrix of indices}

\item{y}{A vector of the non-bootstrapped outcomes}

\item{multiplicity}{A vector of the multiplicity for each nearest neighbors observation}

\item{data}{Vector of the bootstrapped indices}

\item{K}{A integer for K}

\item{self_as_NN}{A boolean where to consider a point itself as its own nearest neighbor}
}
\value{
A list of the weighted sum and number of NN
}
\description{
This function computes the cumulative weighted sum over K nearest neighbors
}
\keyword{internal}
