#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat safe_entropy(const arma::mat& X) {
  arma::mat A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}

// [[Rcpp::export]]
double soft_kmeans_score_cpp(const arma::mat& R, const arma::rowvec& w, const arma::mat& dist_mat, float sigma) {
    float score_dist = arma::as_scalar(arma::accu((R * arma::diagmat(w)) % dist_mat)); 
//     float score_dist = arma::as_scalar(arma::accu((R.each_row() % w) % dist_mat)); 
    float score_entropy = arma::as_scalar(arma::accu((safe_entropy(R)) * arma::diagmat(w)));
    return score_dist + sigma*score_entropy;
}

bool kmeans_converged(const vector<float>& scores, float tol) {
    float s0 = scores[scores.size() - 2];
    float s1 = scores[scores.size() - 1];
    return (s0 - s1) / s0 < tol;
}

// [[Rcpp::export]]
List soft_kmeans_weighted_cpp(arma::mat Y, arma::mat Z, const arma::rowvec& w, unsigned max_iter, float sigma, float tol) {
    Y = arma::normalise(Y, 2, 0); // L2 normalize the columns
    Z = arma::normalise(Z, 2, 0); // L2 normalize the columns
    arma::mat R, Rw, dist_mat; 
    std::vector<float> scores;
    float s0, s1;
    for (unsigned i = 0; i < max_iter; i++) {
        dist_mat = 2 * (1 - Y.t() * Z); 
        R = -dist_mat / sigma;
        R.each_row() -= arma::max(R, 0);  
        R = exp(R);
        R.each_row() /= arma::sum(R, 0);
        Rw = R * arma::diagmat(w); 
        scores.push_back(soft_kmeans_score_cpp(R, w, dist_mat, sigma));
        Y = arma::normalise(Z * Rw.t(), 2, 0); 
        if (i > 0 && kmeans_converged(scores, tol)) break;
    }

    List result = List::create(_("R") = R , _["Y"] = Y, _["scores"] = scores);
    return result;    
}



