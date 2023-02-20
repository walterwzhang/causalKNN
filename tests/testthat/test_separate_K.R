context("compare.separate.results")



test_that("Separate KNN Approach Works", {

    ## Generate the Data  -------------------------------------------------------------------------

    set.seed(1234L)
    n     <-   5000  # Number of Observations
    p     <-   21    # Number of Features
    key   <-   1:n   # Individual indices
    X     <-   cbind(matrix(rnorm(p/3*n), n, p/3),
                     matrix(5*runif(p/3*n), n, p/3),
                     matrix(rexp(p/3*n), n, p/3))
    X_test           <-   matrix(0, 501, p)
    X_test[,1]       <-   rnorm(501)
    X_test[,p/3+1]   <-   5*runif(501)
    TE   <-   pmin(abs(X[,1] * X[,p/3+1] + 1), 4)
    W    <-   rbinom(n, 1, 0.5)
    Y    <-   TE * W + 5 * X[,2] + X[,3] * pmax(X[,5], 0) + 5 * sin(X[,4]) + X[,8] + rnorm(n)

    ##  -------------------------------------------------------------------------------------------



    ## Sequential Run -----------------------------------------------------------------------------

    k                 <-   150L         # Number of NN to search
    distance_metric   <-   "euclidean"  # Distance metric
    keep_dist         <-   FALSE        # Drop the distance matricies
    threads           <-   1L           # Number of threads to use

    knn_index_list    <-   knn_index_mat(X, W, k,
                                        distance_metric = distance_metric,
                                        keep_dist       = keep_dist,
                                        threads         = threads)

    N_step   <-   25  # How many steps to take
    K_step   <-   5   # Value between steps and starting value

    optimal_k_list   <-   knn_optimal_k_separate(X, W, Y, key, N_step, K_step, knn_index_list)

    optimal_K_values   <-   c(optimal_k_list$optimal_K_untreated, optimal_k_list$optimal_K_treated)
    causal_KNN_DF      <-   causalknn_treatment_effect(W, Y, key,
                                                       optimal_K_values, knn_index_list)

    TEP_results_list   <-   tep_enet(X, causal_KNN_DF$treatment_effect, X_test, alpha = 0)

    ##  -------------------------------------------------------------------------------------------

    expect_true(TRUE)
})



test_that("Separate KNN with same K value is reasonable", {

    ## Generate the Data  -------------------------------------------------------------------------

    set.seed(1234L)
    n     <-   500   # Number of Observations
    p     <-   21    # Number of Features
    key   <-   1:n   # Individual indices
    X     <-   cbind(matrix(rnorm(p/3*n), n, p/3),
                     matrix(5*runif(p/3*n), n, p/3),
                     matrix(rexp(p/3*n), n, p/3))
    X_test           <-   matrix(0, 501, p)
    X_test[,1]       <-   rnorm(501)
    X_test[,p/3+1]   <-   5*runif(501)
    TE   <-   pmin(abs(X[,1] * X[,p/3+1] + 1), 4)
    W    <-   rbinom(n, 1, 0.5)
    Y    <-   TE * W + 5 * X[,2] + X[,3] * pmax(X[,5], 0) + 5 * sin(X[,4]) + X[,8] + rnorm(n)

    ##  -------------------------------------------------------------------------------------------



    ## Sequential Run -----------------------------------------------------------------------------


    k                 <-   150L         # Number of NN to search
    distance_metric   <-   "euclidean"  # Distance metric
    keep_dist         <-   FALSE        # Drop the distance matricies
    threads           <-   1L           # Number of threads to use

    knn_index_list    <-   knn_index_mat(X, W, k,
                                        distance_metric = distance_metric,
                                        keep_dist       = keep_dist,
                                        threads         = threads)

    N_step   <-   25  # How many steps to take
    K_step   <-   5   # Value between steps and starting value

    # Suppress max K reached warning - not testing for that here
    suppressWarnings(optimal_k_list  <- knn_optimal_k(X, W, Y, key, N_step, K_step, knn_index_list))

    causal_KNN_DF       <-   causalknn_treatment_effect(W, Y, key,
                                                        optimal_k_list$optimal_K, knn_index_list)

    # Same K for treated and untreated
    optimal_K_values    <-   c(optimal_k_list$optimal_K,
                               optimal_k_list$optimal_K)
    causal_KNN_DF_sep   <-   causalknn_treatment_effect(W, Y, key,
                                                        optimal_K_values, knn_index_list)


    ##  -------------------------------------------------------------------------------------------



    ## Test for Equality --------------------------------------------------------------------------

    expect_equal(causal_KNN_DF, causal_KNN_DF_sep)

    ##  -------------------------------------------------------------------------------------------
})
