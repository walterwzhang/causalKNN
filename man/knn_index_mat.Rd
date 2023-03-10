% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_index_mat.R
\name{knn_index_mat}
\alias{knn_index_mat}
\title{Compute the KNN index matrix}
\usage{
knn_index_mat(
  DF,
  W,
  k = floor(sqrt(nrow(DF))),
  DF_test = NULL,
  distance_metric = "euclidean",
  standardize = TRUE,
  keep_dist = FALSE,
  threads = 1L,
  cov_DF = NULL
)
}
\arguments{
\item{DF}{A data frame of the features in the training sample                             (data.frame)}

\item{W}{A vector of the treatment indicator (1/0 coded) in the training sample           (integer)}

\item{k}{A integer for the max k value for the of kNN. Defaults to floor(sqrt(nrow(DF)))  (integer)}

\item{DF_test}{A data frame of the features in the test sample                            (data.frame)}

\item{distance_metric}{String supplying the method                                        (character)}

\item{standardize}{A boolean flag for whether to standardize the features                 (logical)}

\item{keep_dist}{A boolean flag for keeping generated distance matrices                   (logical)}

\item{threads}{A integer determining the number of threads used - uses openMP             (integer)}

\item{cov_DF}{A custom supplied covariance matrix for the mahalanobis distance metric     (data.frame)}
}
\value{
A list containing the index (and distance) matrices                               (list)
}
\description{
Find the nearest neighbor indices for the treated and untreated for a supplied data.frame
Uses the knn.index.dist function from the KernelKnn package
Currently only applicable to a binary treatment coded as (1/0)
In the returned matrices, the rows are indexed by NN, and the columns are the DF observations
}
