# knn_index_mat.R
# Contains function that find the nearest-neighbor index matrix
#
# Key functions:
#   - knn_index_mat



# Functions =======================================================================================

# knn_index_mat -----------------------------------------------------------------------------------
#' Compute the KNN index matrix
#'
#' Find the nearest neighbor indices for the treated and untreated for a supplied data.frame
#' Uses the knn.index.dist function from the KernelKnn package
#' Currently only applicable to a binary treatment coded as (1/0)
#' In the returned matrices, the rows are indexed by NN, and the columns are the DF observations
#'
#' @param DF A data frame of the features in the training sample                             (data.frame)
#' @param W A vector of the treatment indicator (1/0 coded) in the training sample           (integer)
#' @param k A integer for the max k value for the of kNN. Defaults to floor(sqrt(nrow(DF)))  (integer)
#' @param DF_test A data frame of the features in the test sample                            (data.frame)
#' @param distance_metric String supplying the method                                        (character)
#' @param standardize A boolean flag for whether to standardize the features                 (logical)
#' @param keep_dist A boolean flag for keeping generated distance matrices                   (logical)
#' @param threads A integer determining the number of threads used - uses openMP             (integer)
#' @param cov_DF A custom supplied covariance matrix for the mahalanobis distance metric     (data.frame)
#' @return A list containing the index (and distance) matrices                               (list)
#' @export
#' @importFrom KernelKnn knn.index.dist

knn_index_mat   <-   function(DF, W,
                              k                 =   floor(sqrt(nrow(DF))),
                              DF_test           =   NULL,
                              distance_metric   =   "euclidean",
                              standardize       =   TRUE,
                              keep_dist         =   FALSE,
                              threads           =   1L,
                              cov_DF            =   NULL)
{
    # Checks
    if (length(W) != nrow(DF))
        stop("Length of the W vector is not the same as the nrows of DF")
    if (threads < 0 | threads %% 1 != 0)
        stop("Incorrect threads specification")
    if (!is.null(DF_test))
          if(ncol(DF) != ncol(DF_test))
            stop("Training and test data have a different number of features")

    # Standardize the Features if needed
    if (standardize)
    {
        DF <- apply(DF, 2, standardizeVector)
        if (!is.null(DF_test))
        {
            DF_test <- apply(DF_test, 2, standardizeVector)
        }
    }

    # Subset Data to treated and untreated
    treated_DF     <-   subset(DF, W == 1L)
    untreated_DF   <-   subset(DF, W == 0L)


    if (is.null(DF_test))
    {
        # No test set case
        knn_test_DF   <-   DF

    } else {
        # Test set case
        knn_test_DF   <-   DF_test
    }

    # Find the Treated and Untreated NN
    if (distance_metric == "mahalanobis" & !is.null(cov_DF))
    {
        # Case with the mahalanobis distance metric with a user supplied covariance matrix
        treated_result     <-   knn_index_dist_cov_mahalanobis(TRAIN_DATA   =   treated_DF,
                                                               TEST_DATA    =   knn_test_DF,
                                                               k            =   k,
                                                               threads      =   threads,
                                                               COV_MAT      =   cov_DF)
        untreated_result   <-   knn_index_dist_cov_mahalanobis(TRAIN_DATA   =   untreated_DF,
                                                               TEST_DATA    =   knn_test_DF,
                                                               k            =   k,
                                                               threads      =   threads,
                                                               COV_MAT      =   cov_DF)
    } else
    {
        # Non-mahalanobis distance metric case
        treated_result     <-   KernelKnn::knn.index.dist(data        =   treated_DF,
                                                          TEST_data   =   knn_test_DF,
                                                          k           =   k,
                                                          method      =   distance_metric,
                                                          threads     =   threads)

        untreated_result   <-   KernelKnn::knn.index.dist(data        =   untreated_DF,
                                                          TEST_data   =   knn_test_DF,
                                                          k           =   k,
                                                          method      =   distance_metric,
                                                          threads     =   threads)
    }

    if (!keep_dist) treated_result$test_knn_dist <- NULL
    treated_result$test_knn_ind   <-   apply(treated_result$test_knn_idx, 2, as.integer)
    if (!keep_dist) untreated_result$test_knn_dist <- NULL
    untreated_result$test_knn_ind   <-   apply(untreated_result$test_knn_idx, 2, as.integer)

    # Transpose the matrices to be rows as NN, and cols as obs
    treated_result    <-  lapply(treated_result, t)
    untreated_result  <-  lapply(untreated_result, t)

    # Re-index the values to reflect the training DF
    index       <-    1:nrow(DF)
    index_0     <-    index[W == 0L]
    index_1     <-    index[W == 1L]

    untreated_result$test_knn_ind   <-   t(matrix(index_0[c(untreated_result$test_knn_ind)],
                                                  ncol = k))
    treated_result$test_knn_ind     <-   t(matrix(index_1[c(treated_result$test_knn_ind)],
                                                  ncol = k))

    # Return a list of the matrices
    if (keep_dist)
    {
        return(list(knn_index_matricies = list(untreated_index_mat   =   untreated_result$test_knn_ind,
                                               treated_index_mat     =   treated_result$test_knn_ind),
                    knn_dist_matricies  = list(untreated_index_mat   =   untreated_result$test_knn_dist,
                                               treated_index_mat     =   treated_result$test_knn_dist)
                ))
    }
    return(list(knn_index_matricies = list(untreated_index_mat   =   untreated_result$test_knn_ind,
                                           treated_index_mat     =   treated_result$test_knn_ind)))
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================



# Helper Functions ================================================================================


# standardizeVector ------------------------------------------------------------------------------
#' Standardize a vector to mean 0 and variance 1
#'
#' @param x A vector to be standardized (vector)
#' @return A standardized vector (vector)
#' @keywords internal
#' @importFrom stats sd
#' @noRd

standardizeVector <- function(x)
{
    mean_x   <-   mean(x, na.rm = TRUE)
    sd_x     <-   sd(x,   na.rm = TRUE)
    if (sd_x != 0)
    {
        return((x - mean_x)/sd_x)
    } else
    {
        return((x - mean_x))
    }

}

# -------------------------------------------------------------------------------------------------



# knn_index_dist_cov_mahalanobis ---------------------------------------------------------------------
#' Wraps the \code{knn_index_dist_rcpp} function.
#'
#' General change from the original package's file is that the mahalanobis distance now allows
#' for a custom co-variance matrix or substitute to be used.
#' Still supports many other distance metrics such as euclidian or chebyshev. See the original
#' "KernelKnn" R Package for more details.
#'
#' Returns a list of NN indicies (matrix) and distances (matrix). For both matricies, the columns
#' equal the NN and the Rows are the observation.
#'
#' @param TRAIN_DATA Matrix of the "training" data (matrix)
#' @param k Number of NN to find (integer)
#' @param TEST_DATA Matrix of the "test" data (matrix)
#' @param threads Number of threads to use (integer)
#' @param COV_MAT User provided covariance matric for mahalanobis distance (matrix)
#' @return A list containing matricies of the nearest neighbors indicies and distances (list)
#' @keywords internal
#' @export

knn_index_dist_cov_mahalanobis   <-   function(TRAIN_DATA, k,
                                               TEST_DATA   =   NULL,
                                               threads     =   1L,
                                               COV_MAT     =   NULL)
{
    stopifnot(class(k) %in% c("numeric", "integer"))
    stopifnot(class(threads) %in% c("numeric", "integer"))
    stopifnot(!is.null(TRAIN_DATA))

    if (!inherits(TRAIN_DATA, "matrix"))
    {
        TRAIN_DATA   <-   as.matrix(TRAIN_DATA)
    }

    if (is.null(TEST_DATA))
    {
        TEST_DATA   <-    matrix(nrow = 0, ncol = 0)
    } else if (!inherits(TEST_DATA, "matrix"))
    {
        TEST_DATA   <-   as.matrix(TEST_DATA)
    }

    if (is.null(COV_MAT))
    {
        COV_MAT   <-    matrix(nrow = 0, ncol = 0)
    } else if (!inherits(COV_MAT, "matrix"))
    {
        COV_MAT   <-   as.matrix(COV_MAT)
    }

    if (inherits(threads, "numeric"))
    {
        threads   <-  as.integer(threads)
    }

    if (inherits(k, "numeric"))
    {
        k   <-  as.integer(k)
    }

    return(knn_index_dist_rcpp(TRAIN_DATA, TEST_DATA, k, threads, COV_MAT))
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================
