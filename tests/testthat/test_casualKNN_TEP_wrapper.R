context("compare.wrapper.results")



test_that("causalKNN_TEP produces same results a sequential run", {

    ## Generate the Data  -------------------------------------------------------------------------

     set.seed(1234L)
    n     <-   5000L  # Number of Observations
    p     <-   21L    # Number of Features
    key   <-   1L:n   # Individual indices
    X     <-   cbind(matrix(rnorm(p/3*n), n, p/3),
                     matrix(5*runif(p/3*n), n, p/3),
                     matrix(rexp(p/3*n), n, p/3))
    n_test           <-   501L
    X_test           <-   matrix(0, n_test, p)
    X_test[,1]       <-   rnorm(n_test)
    X_test[,p/3+1]   <-   5*runif(n_test)
    key_test         <-    1:n_test
    TE   <-   pmin(abs(X[,1] * X[,p/3+1] + 1), 4)
    W    <-   rbinom(n, 1, 0.5)
    Y    <-   TE * W + 5 * X[,2] + X[,3] * pmax(X[,5], 0) + 5 * sin(X[,4]) + X[,8] + rnorm(n)

    ##  -------------------------------------------------------------------------------------------



    ## Sequential Run -----------------------------------------------------------------------------
    k                 <-   120L         # Number of NN to search
    distance_metric   <-   "euclidean"  # Distance metric
    keep_dist         <-   FALSE        # Drop the distance matricies
    threads           <-   1L           # Number of threads to use

    knn_index_list    <-   knn_index_mat(X, W, k,
                                        distance_metric = distance_metric,
                                        keep_dist       = keep_dist,
                                        threads         = threads)

    N_step   <-   20L  # How many steps to take
    K_step   <-   5L   # Value between steps and starting value

    optimal_k_list   <-   knn_optimal_k(X, W, Y, key, N_step, K_step, knn_index_list)

    causal_KNN_DF   <-   causalknn_treatment_effect(W, Y, key, optimal_k_list$optimal_K, knn_index_list)

    TEP_results_list   <-   tep_enet(X, causal_KNN_DF$treatment_effect, X_test, alpha = 0)

    ##  -------------------------------------------------------------------------------------------



    ## Wrapper Function ---------------------------------------------------------------------------

    results_list   <-   causalKNN_TEP(X, W, Y, key, X_test, key_test,
                                                            knn_index_mat_options = list(k = 120L),
                                                            knn_optimal_k_options = list(N_step = 20L,
                                                                                         K_step = 5L),
                                                            tep_e_net_options = list(alpha = 0))

    ##  -------------------------------------------------------------------------------------------



    ## Equality Check -----------------------------------------------------------------------------

    expect_equal(TEP_results_list$test_TE, results_list$TEP_results$test_TE)

    ##  -------------------------------------------------------------------------------------------
})
