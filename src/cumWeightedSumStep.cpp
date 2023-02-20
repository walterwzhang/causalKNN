/*  cumWeightedSumStep.cpp
 *  ------------------------------------------------------------------------------------------------
 *  - Computes the KNN Treatment Effect for a vector of K values (defined by N_steps and K_step)
 */



#include <Rcpp.h>
using namespace Rcpp;



//' cumWeightedSumStep
//'
//' This function computes the cumulative weighted sum over K nearest neighbors for a sequence of K values
//'
//' @param Index The matrix of nearest neighbor indices
//' @param y A vector of the non-bootstrapped outcomes
//' @param multiplicity A vector of the multiplicity for each nearest neighbors observation
//' @param data Vector of the bootstrapped indices
//' @param N_steps A integer for the number of steps of K_step to take
//' @param K_step A integer for the value of steps of K and the initial K value
//' @return A list of the weighted sum and number of NN
//' @export
//' @keywords internal
// [[Rcpp::export]]
List cumWeightedSumStep(IntegerMatrix Index, NumericVector y, IntegerVector multiplicity, IntegerVector data,
                        const int N_steps, const int K_step) {



   const int L = data.size();             // Number of observations in data index vector
   const int D = Index.nrow();            // Number of nearest neighbors in index matrix

   NumericMatrix x(N_steps, L);
   IntegerMatrix cum_n(N_steps, L);

   const int K = N_steps*K_step;
   NumericVector v(K);
   IntegerVector mult(K);


   for (int i = 0; i < L; i++) {
      const int index = data[i] - 1;

      int l        = 0;       // Index in v
      int k        = 0;       // Index in the index matrix
      int cum_mult = 0;

      while (cum_mult < K && k < D) {
         int j = Index(k, index) - 1;

         if (j == index)
         {
            k++;       // Skips Self as NN
            continue;
         }

         int m = multiplicity[j];

         if (cum_mult + m > K) m = K - cum_mult;


         for (int n = 0; n < m; n++) {
            v[l]    = y[j];
            mult[l] = 1;
            l++;
         }

         cum_mult += m;
         k++;
      }

      for (int s = 0; s < N_steps; s++) {

         if (s > 0) {
            x(s, i)     = x(s-1, i);
            cum_n(s, i) = cum_n(s-1, i);
         }

         for (int n = (s*K_step); n < ((s+1)*K_step); n++) {
            x(s, i)    += v[n];
            cum_n(s, i) += mult[n];
         }
      }
   }


   return List::create(Named("weighted_sum") = x,
                       Named("cum_n")        = cum_n);
}
