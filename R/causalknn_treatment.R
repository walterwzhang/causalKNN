# causalknn_treatment.R
# Computes the causal KNN Treatment Effect with a given K and KNN Index Matrices
#
# Key Functions:
#     - causalknn_treatment



# Functions =======================================================================================

# causalknn_treatment -----------------------------------------------------------------------------
#' Treatment effect estimation with causal KNN for K values.
#'
#' In separate \code{K} value case, supply a vector of two elements to the \code{K} parameter. The first
#' element of vector will be \code{K} for untreated nearest neighbors, and the second elements will be
#' the K for treated nearest neighbors. Supply a single value for K otherwise.
#'
#' Supply a \code{key_test} if looking to estimate treatment effects for the validation set. Otherwise,
#' the default value of NULL will instruct the function to estimate treatment effects for the training set.
#'
#' Allows for bootstrapped case - see the bootstrap vignette.
#'
#' @param W A vector of the treatment indicator (1/0 coded) in the training set    (integer)
#' @param Y A vector of the outcome values in the training set                     (numeric)
#' @param key A vector of the indices of the training set                          (integer)
#' @param K Value(s) for K                                                         (integer)
#' @param knn_index_list A list of the KNN index matrices provided by the knn_index_mat function (list)
#' @param key_test A vector of the indices of the test set (Defaults to NULL)      (integer)
#' @param return_outcomes Whether to return predicted potential (treated/untreated) outcomes  (logical)
#' @param bootstrap_keys A vector of the bootstrapped keys for the training set (default to NULL for no bootstrap case)  (integer)
#' @param bootstrap_keys_test A vector of the bootstrapped keys for the test set (default to NULL for no bootstrap case)  (integer)
#' @return A data.frame of the predicted treatment effect, treatment outcome, and untreated outcome by key (data.frame)
#' @export

causalknn_treatment_effect   <-   function(W, Y, key, K, knn_index_list,
                                           key_test = NULL,
                                           return_outcomes = FALSE,
                                           bootstrap_keys = NULL,
                                           bootstrap_keys_test = NULL)
{
    # Checks
    if (length(K) > 2)
        stop("Incorrect specification for K")

    if (length(K) == 1)
    {
        # One K case
        K_0   <-    K  # Untreated K
        K_1   <-    K  # Treated K
    } else
    {
        # Separate K case
        K_0   <-    K[1]  # Untreated K
        K_1   <-    K[2]  # Treated K
    }

    # Estimate the causal KNN Treatment effect
    keys_standard   <-   list(full_keys      =   key,
                              full_outputs   =   Y)

    if(is.null(key_test))
    {
        # Case for getting treatment effect for training sample
        if (is.null(bootstrap_keys))
        {
          # No Bootstrap case
          bootstrap_keys   <-   key
        }

        knn_estimate_0   <-   bootstrapIndexMatrixProcessor(keys_standard,
                                                            knn_index_list$knn_index_matricies$untreated_index_mat,
                                                            bootstrap_keys, K_0)
        knn_estimate_1   <-   bootstrapIndexMatrixProcessor(keys_standard,
                                                            knn_index_list$knn_index_matricies$treated_index_mat,
                                                            bootstrap_keys, K_1)
        key_out   <-   bootstrap_keys
    } else
    {
        # Case for getting treatment effect for test sample
        keys_standard$new_keys   <-   key_test  # Unbootstrapped test set keys

        if (is.null(bootstrap_keys))
        {
          # No Bootstrap training set case
          bootstrap_keys   <-   key
        }
        if (is.null(bootstrap_keys_test))
        {
          # No Bootstrap test set case
          bootstrap_keys_test   <-   key_test
        }

        knn_estimate_0   <-   bootstrapIndexMatrixProcessorTest(keys_standard,
                                                                knn_index_list$knn_index_matricies$untreated_index_mat,
                                                                bootstrap_keys, bootstrap_keys_test, K_0)
        knn_estimate_1   <-   bootstrapIndexMatrixProcessorTest(keys_standard,
                                                                knn_index_list$knn_index_matricies$treated_index_mat,
                                                                bootstrap_keys, bootstrap_keys_test, K_1)
        key_out   <-   bootstrap_keys_test
    }

    if (return_outcomes)
    {
        # Returns the potential outcomes
        kvalue_treatment_subset_DF   <-   data.frame(key                 =   key_out,
                                                     treatment_effect    =   knn_estimate_1 - knn_estimate_0,
                                                     treated_outcome     =   knn_estimate_1,
                                                     untreated_outcome   =   knn_estimate_0)
    } else
    {
        # Returns the treatment effect only
        kvalue_treatment_subset_DF   <-   data.frame(key                 =   key_out,
                                                     treatment_effect    =   knn_estimate_1 - knn_estimate_0)
    }

    return(kvalue_treatment_subset_DF)
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================



# Helper Functions ================================================================================

# bootstrapIndexMatrixProcessor -------------------------------------------------------------------
#' Finds the NN in the training set for a observation in the training set with a KNN index matrix
#'
#' Takes an pre-searched index matrix and finds the correct nearest nearest neighbors for each
#' bootstrap iteration.
#' In the unbootstrapped run, the keys_standard will be the same as the keys.
#' In a bootstrapped run, the keys parameter will contained the bootstrapped keys.
#'
#' @param keys_standard A list of the unbootstrapped outcomes and keys in the training set (list)
#' @param index_matrix  A matrix of indices of the NN (matrix)
#' @param keys A vector of the keys in the training set (integer)
#' @param k_value A value for K for the causal KNN (integer)
#' @param self_as_NN Boolean flag to use yourself as NN. Defaults to FALSE (logical)
#' @return A vector of KNN estimates (numeric)
#' @export
#' @import data.table
#' @keywords internal

bootstrapIndexMatrixProcessor   <-   function(keys_standard, index_matrix,
                                              keys, k_value,
                                              self_as_NN = FALSE)
{
    index      <-   NULL
    ordering   <-   NULL

    original_keys   <-   keys_standard$full_keys
    outputs         <-   keys_standard$full_outputs

    old_map_DT       <-   data.table::data.table(keys = original_keys, index = 1:length(original_keys))
    new_keys_DT      <-   data.table::data.table(keys = keys, ordering = seq(1, length(keys), 1))
    t1               <-   merge(new_keys_DT, old_map_DT, all.x = TRUE, by = "keys")[order(ordering)]
    bootstrap_data   <-   as.integer(t1[, index])

    f_table                                    <-   table(bootstrap_data)
    multiplicity                               <-   rep(0L, length(original_keys))
    multiplicity[as.integer(names(f_table))]   <-   as.integer(f_table)

    cwsum_result   <-   cumWeightedSum(index_matrix, outputs, multiplicity, bootstrap_data, as.integer(k_value), self_as_NN)
    knn_estimate   <-   cwsum_result$weighted_sum/cwsum_result$cum_n

    # Warning for insufficient NN in the KNN
    if (sum(cwsum_result$cum_n < k_value) > 0)
    {
        message   <-   paste("Number of observations with insufficient number of neighbors:",
                             sum(cwsum_result$cum_n < k_value), "\n")
        warning(message, immediate. = TRUE)
    }

    return(knn_estimate)
}

# -------------------------------------------------------------------------------------------------



# bootstrapIndexMatrixProcessorTest ------------------------------------------------------------------
#' Finds the NN in the training set for a observation in the test set with a KNN index matrix
#'
#' Takes an pre-searched index matrix and finds the correct nearest nearest neighbors for each
#' bootstrap iteration.
#' In the unbootstrapped run, the keys_standard will be the same as the keys.
#' In a bootstrapped run, the keys parameter will contained the bootstrapped keys
#'
#' @param keys_standard A list of the unbootstrapped outcomes (training set) and keys (training and test sets) (list)
#' @param index_matrix  A matrix of indices of the NN (matrix)
#' @param keys A vector of the keys in the training set (integer)
#' @param keys_test A vector of the keys in the test set (integer)
#' @param k_value A value for K for the causal KNN (integer)
#' @return A vector of KNN estimates (numeric)
#' @export
#' @import data.table
#' @keywords internal

bootstrapIndexMatrixProcessorTest   <-   function(keys_standard, index_matrix,
                                                  keys, keys_test, k_value)
{
  index      <-   NULL
  ordering   <-   NULL

  original_keys   <-   keys_standard$full_keys
  outputs         <-   keys_standard$full_outputs
  new_keys        <-   keys_standard$new_keys

  old_map_DT       <-   data.table::data.table(keys = original_keys, index = 1:length(original_keys))
  new_keys_DT      <-   data.table::data.table(keys = keys, ordering = seq(1, length(keys), 1))
  t1               <-   merge(new_keys_DT, old_map_DT, all.x = TRUE, by = "keys")[order(ordering)]
  bootstrap_data   <-   as.integer(t1[, index])

  f_table                                    <-   table(bootstrap_data)
  multiplicity                               <-   rep(0L, length(original_keys))
  multiplicity[as.integer(names(f_table))]   <-   as.integer(f_table)

  new_map_DT             <-   data.table::data.table(keys = new_keys, index= 1:length(new_keys))
  new_keys_valid_DT      <-   data.table::data.table(keys = keys_test, ordering = seq(1, length(keys_test), 1))
  t1_valid               <-   merge(new_keys_valid_DT, new_map_DT, all.x = TRUE, by = "keys")[order(ordering)]
  bootstrap_valid_data   <-   as.integer(t1_valid[, index])

  cwsum_result   <-   cumWeightedSum(index_matrix, outputs, multiplicity, bootstrap_valid_data, as.integer(k_value))
  knn_estimate   <-   cwsum_result$weighted_sum/cwsum_result$cum_n

  # Warning for insufficient NN in the KNN
  if (sum(cwsum_result$cum_n < k_value) > 0)
  {
      message   <-   paste("Number of observations with insufficient number of neighbors:",
                           sum(cwsum_result$cum_n < k_value), "\n")
      warning(message, immediate. = TRUE)
  }

  return(knn_estimate)
}

# -------------------------------------------------------------------------------------------------

# =================================================================================================
