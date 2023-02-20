# causalknn_tep.R
# Wraps the knn_index_mat, knn_optimal_k, causalknn_treatment_effect, and tep_projection_e_net
# functions to provide an "off-the-shelf" implementation for the Causal KNN TEP.
#
# Key Functions:
#     - causalKNN_TEP




# Functions ========================================================================================

# causalKNN_TEP ------------------------------------------------------------------------------------
#' Performs Causal KNN Regression and Treatment Effect Projection
#'
#' Wraps the causal KNN Regression and the Treatment Effect Projection into one function. The
#' treatment effect projection step is performed with a Elastic-Net
#' Sequentially runs \code{knn_index_mat}, \code{knn_optimal_k},
#' \code{causalknn_treatment_effect},  and \code{tep_projection_e_net} functions and returns a list
#' of treatment results and optimal parameters. Since all of the important parameters are returned,
#' the user can re-run certain functions easily. (i.e. re-running \code{causalknn_treatment_effect} with
#' a larger K value).
#'
#' This off-the-shelf implementation does not implement bootstrapping for finding the optimal K
#' value and for the treatment effect projection.
#'
#' @param DF A data frame of the features                                                                        (data.frame)
#' @param W A vector of the treatment indicator (1/0 coded)                                                      (integer)
#' @param Y A vector of the outcome values                                                                       (numeric)
#' @param key A vector of the indices                                                                            (integer)
#' @param DF_test A data frame of the features in the test set                                                   (data.frame)
#' @param key_test A vector of the indices of the test sample                                                    (integer)
#' @param threads The value of number of threads to use                                                          (integer)
#' @param knn_index_mat_options A list of parameters to pass into \code{knn_index_mat}                           (list)
#' @param knn_optimal_k_options A list of parameters to pass into \code{knn_optimal_k}                           (list)
#' @param tep_e_net_options A list of parameters to pass into \code{tep_projection_e_net}             (list)
#' @param verbose A boolean flag for verbose output                                                              (logical)
#' @return A list containing function results (list)
#' @export

causalKNN_TEP   <-   function(DF, W, Y, key, DF_test, key_test,
                              threads                              =   1L,
                              knn_index_mat_options                =   list(k = floor(sqrt(nrow(DF)))),
                              knn_optimal_k_options                =   list(N_step = ifelse(knn_index_mat_options$k > 300L,
                                                                                            floor(knn_index_mat_options$k/25),
                                                                                            14L),
                                                                            K_step = ifelse(knn_index_mat_options$k > 300L,
                                                                                            25L,
                                                                                            floor(knn_index_mat_options$k/14))),
                              tep_e_net_options         =   list(),
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

    # Generate KNN Index Matrix
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

    # Causal KNN TE estimate for a specific K value
    causal_KNN_DF   <-   do.call(causalknn_treatment_effect, c(list(W                =   W,
                                                                    Y                =   Y,
                                                                    key              =   key,
                                                                    K                =   optimal_K,
                                                                    knn_index_list   =   knn_index_list)))
    if (verbose) print(paste0("Causal KNN Treatment Effect Estimated. Elapsed Time: ", format(Sys.time() - init_time)))

    # TEP Projection with the Elastic Net
    TEP_results   <-   do.call(tep_enet, c(list(training_X    =   DF,
                                                training_TE   =   causal_KNN_DF$treatment_effect,
                                                test_X        =   DF_test),
                                           tep_e_net_options))

    if (verbose) print(paste0("Treatment Effect Projection Complete.  Elapsed Time: ", format(Sys.time() - init_time)))

    # Create a data frame of the estmated treatment effects
    causal_KNN_test_DF   <-   data.frame(key                =   key_test,
                                         treatment_effect   =   TEP_results$test_TE)

    # Return the results from the four functions in a list
    return(list(knn_index_list        =   knn_index_list,
                optimal_k_list        =   optimal_k_list,
                causal_KNN_DF         =   causal_KNN_DF,
                TEP_results           =   TEP_results,
                causal_KNN_test_DF    =   causal_KNN_test_DF))
}

# --------------------------------------------------------------------------------------------------

# ==================================================================================================
