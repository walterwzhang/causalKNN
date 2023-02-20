/*  mahalanobis_covariance_distance.cpp
 *  ------------------------------------------------------------------------------------------------
 *
 *  Key Functions:
 *    - knn_index_dist_rcpp
 */



#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif



// Helper Functions ===============================================================================

// INV_EXC ----------------------------------------------------------------------------------------
//' INV_EXC
//'
//' Computes the inverse or the Moore-Penrose pseudo-inverse when \code{supplied} is true. If
//' \code{supplied} is false, then it computes the covariance matrix and then it's inverse.
//'
//' @param cov_data A matrix of the data or the covariance matrix
//' @param supplied A boolean of whether the covariance matrix is supplied.
//' @return The inverse covariance matrix
//' @noRd
// [[Rcpp::export]]
arma::mat INV_EXC(arma::mat cov_data, bool supplied = false)
{

    arma::mat inverse_temp;

    try {
        if (supplied) {
            inverse_temp = arma::inv(cov_data);
        } else {
            inverse_temp = arma::inv(arma::cov(cov_data));
        }
    } catch (...) {
        Rcpp::warning("Pseudo-inverse is calculated due to singular input matrix.");
    }

    if (inverse_temp.empty()) {
        if (supplied) {
            inverse_temp = arma::pinv(cov_data);
        } else
        {
            inverse_temp = arma::pinv(arma::cov(cov_data));
        }
    }

    return inverse_temp;
}

// ------------------------------------------------------------------------------------------------




// inner_loop -------------------------------------------------------------------------------------
//' inner_loop
//'
//' Computes the inner loop for the Nearest Neighbor computations
//'
//' @param matrix_1 A matrix of data
//' @param matrix_2 A matrix of data
//' @param i A integer denoting the row index
//' @param row_num A integer denoting the number of rows of matrix_2
//' @param threads A integer for the number of threads to use
//' @param cov_mat A supplied covariance matrix
//' @return A row vector of the computed differences with the Mahalanobis distance metric
//' @noRd
// [[Rcpp::export]]
arma::rowvec inner_loop(arma::mat& matrix_1, arma::mat& matrix_2,
                        int i, int row_num, int threads,
                        arma::mat cov_mat)
{

    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #endif

    double temp_ind;

    arma::rowvec temp_out = arma::zeros<arma::rowvec>(row_num);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif

    for (int j = 0; j < row_num; j++) {

        // Only for method == "mahalanobis" and takes in a pre-created covariance matrix
        temp_ind = arma::as_scalar(
                        std::sqrt(arma::as_scalar(
                            ((matrix_1.row(j) - matrix_2.row(i)) * cov_mat) *
                             (matrix_1.row(j) - matrix_2.row(i)).t()
                            )));

        if (temp_ind != temp_ind) {

            // distance if 1.0 for NAs
            // Should not occur if the data is pre-processed correctly
            temp_out(j) = 1.0;
        } else {

            temp_out(j) = temp_ind;
        }
    }

    return temp_out;
}

// ------------------------------------------------------------------------------------------------



// train_mat --------------------------------------------------------------------------------------
//' train_mat
//'
//' Finds the k Nearest Neighbors for the matrix
//'
//' @param MATRIX A matrix of the data
//' @param k An integer of the nearest neighbors to find
//' @param threads An integer of threads to use for Open MP
//' @param COV_MAT The supplied covariance matrix
//' @return A list of the k nearest neighbor indices and distances
//' @noRd
// [[Rcpp::export]]
Rcpp::List train_mat(arma::mat& MATRIX, int k,
                       int threads, arma::mat& COV_MAT)
{

    arma::mat knn_indices;
    arma::mat knn_distances;
    arma::mat cov_mat;

    // For method == "mahalanobis"
    if (COV_MAT.is_empty())
    {
        cov_mat = INV_EXC(arma::cov(MATRIX));
    } else {
        cov_mat = INV_EXC(COV_MAT, true);
    }

    int ITERS = MATRIX.n_rows;

    knn_indices.set_size(ITERS, k);
    knn_distances.set_size(ITERS, k);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif

    for (int i = 0; i < ITERS; i++) {

        arma::rowvec temp_out    =   inner_loop(MATRIX, MATRIX, i, ITERS, threads, cov_mat);
        arma::rowvec index_out   =   arma::conv_to<arma::rowvec>::from(arma::sort_index(temp_out, "ascend")) + 1;
        knn_indices.row(i)       =   index_out.subvec(1, k);
        arma::rowvec temp_sort   =   arma::conv_to<arma::rowvec>::from(arma::sort(temp_out, "ascend"));
        knn_distances.row(i)     =   temp_sort.subvec(1, k);
    }

    return Rcpp::List::create(knn_indices, knn_distances);

}

// ------------------------------------------------------------------------------------------------



// test_mat ---------------------------------------------------------------------------------------
//' test_mat
//'
//' Finds k Nearest Neighbors of the Test Data in the Training Data
//'
//' @param MATRIX A matrix of the training data
//' @param TEST_DATA A matrix of the test data
//' @param k An integer of the nearest neighbors to find
//' @param threads An integer of threads to use for Open MP
//' @param COV_MAT The supplied covariance matrix
//' @return A list of the k nearest neighbor indices and distances
//' @noRd
// [[Rcpp::export]]
Rcpp::List test_mat(arma::mat& MATRIX, arma::mat& TEST_DATA,
               int k, int threads, arma::mat& COV_MAT)
{

    arma::mat knn_indices;
    arma::mat knn_distances;
    arma::mat cov_mat;

    // For method == "mahalanobis"
    if (COV_MAT.is_empty())
    {
        cov_mat = INV_EXC(arma::cov(arma::join_vert(MATRIX, TEST_DATA)));
    } else {

        cov_mat = INV_EXC(COV_MAT, true);
    }


    int ITER_TEST  = TEST_DATA.n_rows;
    int ITER_TRAIN = MATRIX.n_rows;
    knn_indices.set_size(ITER_TEST, k);
    knn_distances.set_size(ITER_TEST, k);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif

    for (int i = 0; i < ITER_TEST; i++) {

        arma::rowvec temp_out    =   inner_loop(MATRIX, TEST_DATA,
                                                i, ITER_TRAIN,
                                                threads,
                                                cov_mat);
        arma::rowvec index_out   =   arma::conv_to<arma::rowvec>::from(arma::sort_index(temp_out, "ascend")) + 1;
        knn_indices.row(i)       =   index_out.subvec(0, k - 1);
        arma::rowvec temp_sort   =   arma::conv_to<arma::rowvec>::from(arma::sort(temp_out, "ascend"));
        knn_distances.row(i)     =   temp_sort.subvec(0, k - 1);
    }

    // return { knn_indices, knn_distances };
    return Rcpp::List::create(knn_indices, knn_distances);
}

// ------------------------------------------------------------------------------------------------

// ================================================================================================



// Main Functions =================================================================================

// knn_index_dist_rcpp ----------------------------------------------------------------------------
//' knn_index_dist_rcpp
//'
//' Compute NN for Mahalanobis Distance Matricies and wraps the Rcpp helper functions
//'
//' @param MATRIX A matrix of training data
//' @param TEST_DATA A matrix of test data
//' @param k A integer of the Nearest Neighbors to find
//' @param threads A integer of the number of threads to use with Open MP
//' @param COV_MAT User supplied covariance matrix
//' @return A list of the computed nearest neighbor index matrix and distance matrix
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List knn_index_dist_rcpp(arma::mat& MATRIX, arma::mat& TEST_DATA,
                               int k, int threads, arma::mat& COV_MAT)
{

    // For method == "mahalanobis" only (implicitly)
    std::string name_ind;
    std::string name_dist;
    Rcpp::List data_out;


    if (TEST_DATA.is_empty()) {

        // Case to find the NN in the training sample only
        data_out    =   train_mat(MATRIX, k, threads, COV_MAT);
        name_ind    =   "train_knn_idx";
        name_dist   =   "train_knn_dist";
    }

    if (!TEST_DATA.is_empty()) {

        // Case to find the NN in the training sample from the test sample
        data_out    =   test_mat(MATRIX, TEST_DATA, k, threads, COV_MAT);
        name_ind    =   "test_knn_idx";
        name_dist   =   "test_knn_dist";
    }

    return Rcpp::List::create(Rcpp::Named(name_ind)  = data_out[0],
                              Rcpp::Named(name_dist) = data_out[1]);
}

// ------------------------------------------------------------------------------------------------

// ================================================================================================

