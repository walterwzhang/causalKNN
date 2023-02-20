# knn_optimal_k.R
# Given a knn nearest neighbor indexed matrix and the data, we find the optimal K value
#
# Key Functions:
#   - knn_optimal_k
#   - knn_optimal_k_separate


# Functions =======================================================================================

# knn_optimal_k -----------------------------------------------------------------------------------
#' Find the optimal K for causal KNN regression
#'
#' Finds the optimal K by choosing a grid of possible K values and finding the K values that has the
#' smallest RMSE to the Transformed Outcome.
#' We can also Bootstrap the procedure to find a more robust optimal K value.
#'
#' @param DF A data frame of the features                                          (data.frame)
#' @param W A vector of the treatment indicator (1/0 coded)                        (integer)
#' @param Y A vector of the outcome values                                         (numeric)
#' @param key A vector of the indices                                              (integer)
#' @param N_step The number of steps to take                                       (integer)
#' @param K_step The value between steps of K and the initial K value              (integer)
#' @param knn_index_list A list of the KNN index matrices provided by the knn_index_mat function  (list)
#' @param propensity A vector of the propensity scores (defaults to RCT setting)   (numeric)
#' @param bootstrap_keys A vector of the bootstrapped keys (default to NULL for no bootstrap case)  (integer)
#' @return A list containing the optimal K value (list)
#' @import data.table
#' @export

knn_optimal_k   <-   function(DF, W, Y, key,
                              N_step, K_step, knn_index_list,
                              propensity = mean(W),
                              bootstrap_keys = NULL)
{
    # Checks
    if (length(W) != nrow(DF))
        stop("Length of the W vector is not the same as the nrows of DF")
    if (length(Y) != nrow(DF))
        stop("Length of the Y vector is not the same as the nrows of DF")
    if (length(key) != nrow(DF))
        stop("Length of the key vector is not the same as the nrows of DF")

    if (is.null(bootstrap_keys))
    {
        # No Bootstrap case
      bootstrap_keys   <-   key
    }

    # Estimate the Treatment Effect for each NN
    keys_standard   <-   list(full_keys      =   key,  # Original keys here
                              full_outputs   =   Y)
    knn_prediction_0   <-   bootstrapKNNStepProcessor(keys_standard,
                                                      knn_index_list$knn_index_matricies$untreated_index_mat,
                                                      bootstrap_keys,
                                                      N_step,
                                                      K_step)
    knn_prediction_1   <-   bootstrapKNNStepProcessor(keys_standard,
                                                      knn_index_list$knn_index_matricies$treated_index_mat,
                                                      bootstrap_keys,
                                                      N_step,
                                                      K_step)
    knn_treatment   <-   data.table::data.table(knn_prediction_1 - knn_prediction_0)

    # Compute the Transformed Outcome
    transformed_outcome   <-   Y * (W - propensity) / (propensity * (1 - propensity))

    # Compute the RMSE for at each K-Value
    rmse_results   <-   unlist(lapply(names(knn_treatment), function(x)
    {
        mean((knn_treatment[, get(x)] - transformed_outcome)^2)
    }))
    names(rmse_results) <- names(knn_treatment)

    # Find the optimal K (yields smallest RMSE)
    optimal_rmse   <-   rmse_results[which.min(rmse_results)[1]]
    optimal_K      <-   as.integer(gsub("KNN_", "", (names(which.min(rmse_results)[1]))))

    if (optimal_K == N_step * K_step)
      warning("Optimal K is the maximum K specified. Consider running the procedure again with a higher maximum possible value of K", immediate. = TRUE)

    return(list(optimal_K      =   optimal_K,
                optimal_rmse   =   optimal_rmse,
                k_vec          =   seq(K_step, N_step * K_step, by = K_step),
                rmse_vec       =   rmse_results))

}

# -------------------------------------------------------------------------------------------------



# knn_optimal_k_separate --------------------------------------------------------------------------
#' Find the optimal K separately for treated/untreated for causal KNN regression
#'
#' Finds the optimal K separately for treated and untreated nearest neighbors by choosing a grid
#' of possible K values and finding the K values that has the smallest RMSE to the Transformed
#' Outcome. The combination of the optimal K value for the treated and the untreated that yields
#' the smallest transformed error loss will be returned.
#'
#' We can also bootstrap the procedure to find a more robust optimal K value.
#'
#' @param DF A data frame of the features                                                           (data.frame)
#' @param W A vector of the treatment indicator (1/0 coded)                                         (integer)
#' @param Y A vector of the outcome values                                                          (numeric)
#' @param key A vector of the indices                                                               (integer)
#' @param N_step The number of steps to take                                                        (integer)
#' @param K_step The value between steps of K and the initial K value                               (integer)
#' @param knn_index_list A list of the KNN index matrices provided by the knn_index_mat function    (list)
#' @param propensity A vector of the propensity scores (defaults to RCT setting)                    (numeric)
#' @param bootstrap_keys A vector of the bootstrapped keys (default to NULL for no bootstrap case)  (integer)
#' @param threads The value of number of threads to use                                             (integer)
#' @return A list containing the optimal K values (list)
#' @import data.table parallel
#' @export

knn_optimal_k_separate   <-   function(DF, W, Y, key,
                                       N_step, K_step, knn_index_list,
                                       propensity = mean(W),
                                       bootstrap_keys = NULL,
                                       threads = 1L)
{
    # Checks
    if (length(W) != nrow(DF))
        stop("Length of the W vector is not the same as the nrows of DF")
    if (length(Y) != nrow(DF))
        stop("Length of the Y vector is not the same as the nrows of DF")
    if (threads < 0 | threads %% 1 != 0)
        stop("Incorrect threads specification")

    if (is.null(bootstrap_keys))
    {
        # No Bootstrap case
        bootstrap_keys   <-   key
    }

    # Estimate the Treatment Effect for each NN
    keys_standard   <-   list(full_keys      =   key,  # Original keys here
                              full_outputs   =   Y)
    knn_prediction_0   <-   bootstrapKNNStepProcessor(keys_standard,
                                                      knn_index_list$knn_index_matricies$untreated_index_mat,
                                                      bootstrap_keys,
                                                      N_step,
                                                      K_step)
    knn_prediction_1   <-   bootstrapKNNStepProcessor(keys_standard,
                                                      knn_index_list$knn_index_matricies$treated_index_mat,
                                                      bootstrap_keys,
                                                      N_step,
                                                      K_step)
    knn_treatment   <-   data.table::data.table(knn_prediction_1 - knn_prediction_0)

    # Compute the Transformed Outcome
    transformed_outcome   <-   Y * (W - propensity) / (propensity * (1 - propensity))

    # Generate the Grid of possible combinations of K values for treated/untreated NN
    K_seq    <-   seq(K_step, K_step * N_step, K_step)
    K_grid   <-   data.table::data.table(expand.grid(K_seq, K_seq))
    setnames(K_grid, c("Untreated", "Treated"))

    # Compute the RMSE for each set of K values
    rmse   <-   function(x, y) {sqrt(mean((x-y)^2))}
    K_0    <-   NULL
    K_1    <-   NULL

    if (threads > 1)
    {
        cl   <-   parallel::makeCluster(threads)
        invisible(parallel::clusterEvalQ(cl, {require(data.table)}))
        parallel::clusterExport(cl, c("K_grid", "knn_prediction_0", "knn_prediction_1", "rmse", "transformed_outcome"))

        rmse_list   <-   parallel::parLapply(cl, 1:nrow(K_grid), function(index)
        {
            K_0   <-   K_grid[index, get("Untreated")] # Untreated K-Value
            K_1   <-   K_grid[index, get("Treated")]   # Treated K-Value

            K0_pred   <-   knn_prediction_0[, paste0("KNN_", K_0)]
            K1_pred   <-   knn_prediction_1[, paste0("KNN_", K_1)]

            knn_TE   <-   K1_pred - K0_pred

            data.table(K_0 = K_0,
                       K_1 = K_1,
                       rmse = rmse(knn_TE, transformed_outcome))
        })
        parallel::stopCluster(cl)
    } else
    {
        rmse_list   <-   lapply(1:nrow(K_grid), function(index)
        {
            K_0   <-   K_grid[index, get("Untreated")] # Untreated K-Value
            K_1   <-   K_grid[index, get("Treated")]   # Treated K-Value

            K0_pred   <-   knn_prediction_0[, paste0("KNN_", K_0)]
            K1_pred   <-   knn_prediction_1[, paste0("KNN_", K_1)]

            knn_TE   <-   K1_pred - K0_pred

            data.table(K_0 = K_0,
                       K_1 = K_1,
                       rmse = rmse(knn_TE, transformed_outcome))
        })
    }

    rmse_DT   <-   data.table::rbindlist(rmse_list)

    optimal_K_untreated   <-   rmse_DT[which.min(rmse_DT$rmse),K_0]
    optimal_K_treated     <-   rmse_DT[which.min(rmse_DT$rmse),K_1]

    if (optimal_K_untreated == N_step * K_step)
        warning("Optimal K for untreated NN is the maximum K specified. Consider running the procedure again with a higher maximum possible value of K", immediate. = TRUE)

    if (optimal_K_untreated == N_step * K_step)
        warning("Optimal K for treated NN is the maximum K specified. Consider running the procedure again with a higher maximum possible value of K", immediate. = TRUE)


    return(list(optimal_K_untreated   =   optimal_K_untreated,
                optimal_K_treated     =   optimal_K_treated,
                optimal_rmse          =   rmse_DT[which.min(rmse_DT$rmse),rmse],
                rmse_DT               =   rmse_DT))
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================



# Helper Functions ================================================================================


# bootstrapKNNStepProcessor -----------------------------------------------------------------------
#' Finds the causal KNN estimates for a specified sequence of K-values
#'
#' Helper function for knn_optimal_k.
#' Takes an pre-searched index matrix and finds the correct nearest nearest neighbors for a sequence
#' of K's.
#' Does not include an observation as it own NN
#' Can be used for bootstrapped data - see the bootstrap vignette
#'
#' @param keys_standard A list of the un-bootstrapped outcomes and keys (list)
#' @param index_matrix  A matrix of indices of the NN (matrix)
#' @param keys A vector of bootstrapped keys (integer)
#' @param N_step The number of steps to take (integer)
#' @param k_step The value between steps of K and initial K value (integer)
#' @return A matrix of the KNN estimates where rows in NN and cols (matrix)
#' @export
#' @import data.table
#' @keywords internal

bootstrapKNNStepProcessor   <-   function(keys_standard, index_matrix, keys, N_step, K_step)
{
    index      <-   NULL
    ordering   <-   NULL

    original_keys   <-   keys_standard$full_keys
    outputs         <-   keys_standard$full_outputs

    old_map_DT       <-   data.table::data.table(keys    =   original_keys,
                                                 index   =   1:length(original_keys))
    new_keys_DT      <-   data.table::data.table(keys       =   keys,
                                                 ordering   =   seq(1, length(keys), 1))
    t1               <-   merge(new_keys_DT, old_map_DT, all.x = TRUE, by = "keys")[order(ordering)]
    bootstrap_data   <-   as.integer(t1[, index])

    f_table                                    <-   table(bootstrap_data)
    multiplicity                               <-   rep(0L, length(original_keys))
    multiplicity[as.integer(names(f_table))]   <-   as.integer(f_table)

    cwsum_result   <-   cumWeightedSumStep(index_matrix, outputs, multiplicity, bootstrap_data, N_step, K_step)

    knn_estimate_list   <-   lapply(1:N_step, function(x)
    {
        cwsum_result$weighted_sum[x,] / cwsum_result$cum_n[x,]
    })

    knn_estimate_matrix   <-   matrix(unlist(knn_estimate_list),
                                      nrow    =   length(bootstrap_data),
                                      byrow   =   FALSE)
    colnames(knn_estimate_matrix)   <-   paste0("KNN_", (1:N_step)*K_step)


    return(knn_estimate_matrix)
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================
