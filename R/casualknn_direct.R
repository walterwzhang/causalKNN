# causalknn_direct.R
# Wraps the knn_index_mat, knn_optimal_k, and causalknn_treatment_effect
# functions to provide an "off-the-shelf" implementation of the direct Causal KNN.
#
# Key Functions:
#     - causalKNN_direct




# Functions ========================================================================================

# causalKNN_direct ------------------------------------------------------------------------------------
#' Performs direct Causal KNN regression
#'
#' Wraps the direct Causal KNN regression into one function.
#' The function sequentially runs \code{knn_index_mat}, \code{knn_optimal_k}, and
#' \code{causalknn_treatment_effect} functions and returns a list
#' of treatment results and optimal parameters. Since all of the important parameters are returned,
#' the user can re-run certain functions easily.
#'
#' This off-the-shelf implementation does not implement bootstrapping for finding the optimal K
#' value and for the treatment effect projection.
#'
#' @param DF A data frame of the features of the training sample                                                  (data.frame)
#' @param W A vector of the treatment indicator of the training sample (1/0 coded)                                (integer)
#' @param Y A vector of the outcome values of the training sample                                                 (numeric)
#' @param key A vector of the indices of the training sample                                                      (integer)
#' @param DF_test A data frame of the features in the test set                                                    (data.frame)
#' @param key_test A vector of the indices of the test sample                                                     (integer)
#' @param threads The value of number of threads to use                                                           (integer)
#' @param knn_index_mat_options A list of parameters to pass into \code{knn_index_mat}                            (list)
#' @param knn_optimal_k_options A list of parameters to pass into \code{knn_optimal_k}                            (list)
#' @param causalknn_treatment_effect_options A list of parameters to pass into \code{causalknn_treatment_effect}  (list)
#' @param verbose A boolean flag for verbose output                                                               (logical)
#' @return A list containing function results (list)
#' @export

causalKNN_direct   <-   function(DF, W, Y, key, DF_test, key_test,
                                 threads                              =   1L,
                                 knn_index_mat_options                =   list(k = floor(sqrt(nrow(DF)))),
                                 knn_optimal_k_options                =   list(N_step = ifelse(knn_index_mat_options$k > 300L,
                                                                                               floor(knn_index_mat_options$k/25),
                                                                                               14L),
                                                                               K_step = ifelse(knn_index_mat_options$k > 300L,
                                                                                               25L,
                                                                                               floor(knn_index_mat_options$k/14))),
                                 causalknn_treatment_effect_options   =   list(),
                                 verbose                              =   FALSE
                                 )
{
    # Checks
    if (length(W) != nrow(DF))
        stop("Length of the W vector is not the same as the nrows of DF")
    if (length(Y) != nrow(DF))
        stop("Length of the Y vector is not the same as the nrows of DF")
    if (length(key) != nrow(DF))
        stop("Length of the key vector is not the same as the nrows of DF")

    if (verbose) init_time <- Sys.time()

    # Generate KNN Index Matrix for finding optimal K
    knn_index_list   <-   do.call(knn_index_mat, c(list(DF        =   DF,
                                                        W         =   W,
                                                        threads   =   threads),
                                                   knn_index_mat_options))
    if (verbose) print(paste0("Index Matrix Computed. Elapsed Time: ", format(Sys.time() - init_time)))

    # Find the Optimal K
    optimal_k_list   <-   do.call(knn_optimal_k, c(list(DF               =   DF,
                                                        W                =   W,
                                                        Y                =   Y,
                                                        key              =   key,
                                                        knn_index_list   =   knn_index_list),
                                                   knn_optimal_k_options))
    optimal_K   <-   optimal_k_list$optimal_K
    if (verbose) print(paste0("Optimal K Found.       Elapsed Time: ", format(Sys.time() - init_time)))

    # Generate KNN Index Matrix for NN in the training set for the test set
    knn_index_test_list    <-   do.call(knn_index_mat, c(list(DF        =   DF,
                                                              DF_test   =   DF_test,
                                                              W         =   W,
                                                              threads   =   threads),
                                                         knn_index_mat_options))


    # Direct Causal KNN TE estimate for a specific K value
    causal_KNN_test_DF   <-   do.call(causalknn_treatment_effect, c(list(W                =   W,
                                                                         Y                =   Y,
                                                                         key              =   key,
                                                                         K                =   optimal_K,
                                                                         knn_index_list   =   knn_index_test_list,
                                                                         key_test         =   key_test),
                                                                    causalknn_treatment_effect_options))
    if (verbose) print(paste0("Direct Causal KNN Treatment Effect Estimated. Elapsed Time: ", format(Sys.time() - init_time)))


    # Return the results from the four functions in a list
    return(list(knn_index_list        =   knn_index_list,
                knn_index_test_list   =   knn_index_test_list,
                optimal_k_list        =   optimal_k_list,
                causal_KNN_test_DF    =   causal_KNN_test_DF))
}

# --------------------------------------------------------------------------------------------------

# ==================================================================================================
