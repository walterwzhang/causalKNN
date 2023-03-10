% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{knn_index_dist_rcpp}
\alias{knn_index_dist_rcpp}
\title{knn_index_dist_rcpp}
\usage{
knn_index_dist_rcpp(MATRIX, TEST_DATA, k, threads, COV_MAT)
}
\arguments{
\item{MATRIX}{A matrix of training data}

\item{TEST_DATA}{A matrix of test data}

\item{k}{A integer of the Nearest Neighbors to find}

\item{threads}{A integer of the number of threads to use with Open MP}

\item{COV_MAT}{User supplied covariance matrix}
}
\value{
A list of the computed nearest neighbor index matrix and distance matrix
}
\description{
Compute NN for Mahalanobis Distance Matricies and wraps the Rcpp helper functions
}
\keyword{internal}
