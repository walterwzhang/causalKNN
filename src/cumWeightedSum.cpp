/*  cumWeightedSum.cpp
 *  ------------------------------------------------------------------------------------------------
 *  - Computes the KNN Treatment Effect a specific K value
 */



#include <Rcpp.h>
using namespace Rcpp;

//' cumWeightedSum
//'
//' This function computes the cumulative weighted sum over K nearest neighbors
//'
//' @param Index The matrix of indices
//' @param y A vector of the non-bootstrapped outcomes
//' @param multiplicity A vector of the multiplicity for each nearest neighbors observation
//' @param data Vector of the bootstrapped indices
//' @param K A integer for K
//' @param self_as_NN A boolean where to consider a point itself as its own nearest neighbor
//' @return A list of the weighted sum and number of NN
//' @export
//' @keywords internal
// [[Rcpp::export]]
List cumWeightedSum(IntegerMatrix Index, NumericVector y,
    IntegerVector multiplicity, IntegerVector data,
    const int K, const bool self_as_NN = false) {

   const int L = data.size();             // Number of observations in data index vector
   const int D = Index.nrow();            // Number of nearest neighbors in Index matrix

   NumericVector x(L);
   IntegerVector cum_n(L);


   for (int i = 0; i < L; i++) {

      const int index = data[i] - 1;

      int k        = 0;
      int cum_mult = 0;
      x[i]         = 0.0;
      cum_n[i]     = 0.0;


      while (cum_mult < K && k < D) {

       const int j = Index(k, index) - 1;

       if (!self_as_NN) {
           if (j == index) {
                k++;       // Skips yourself as a Nearest Neighbor
                continue;
            }
        }
        int m = multiplicity[j];

        cum_mult += m;
        if (cum_mult > K) {
            m -= cum_mult - K;
            cum_mult = K;
        }
        x[i] += m*y[j];
        k++;
    }

    cum_n[i] = cum_mult;
}

return List::create(Named("weighted_sum") = x,
                     Named("cum_n")       = cum_n);
}
