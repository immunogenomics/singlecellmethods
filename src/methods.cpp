#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//typedef arma::mat MATTYPE;
//typedef arma::vec VECTYPE;
typedef arma::fmat MATTYPE;
typedef arma::fvec VECTYPE;



// [[Rcpp::export]]
MATTYPE exp_mean(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, int ncol, int nrow, const uvec& groups, const uvec& group_sizes) {
    int ngroups = group_sizes.n_elem;
    MATTYPE res = zeros<MATTYPE>(nrow, ngroups);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            res(i[j], groups[c]) += std::expm1(x[j]);
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) /= group_sizes[c];
        }
    }
        
    return(res);
}



// [[Rcpp::export]]
MATTYPE log_vmr(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, 
                  int ncol, int nrow, const MATTYPE& means,
                  const uvec& groups, const uvec& group_sizes) {
    
    int ngroups = group_sizes.n_elem;
    MATTYPE res = zeros<MATTYPE>(nrow, ngroups);
    MATTYPE nnzero = zeros<MATTYPE>(nrow, ngroups);
    double tmp;
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            tmp = std::expm1(x[j]) - means(i[j], groups(c));
            res(i[j], groups[c]) += tmp * tmp;
            nnzero(i[j], groups(c))++;
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) += (group_sizes[c] - nnzero(r, c)) * means(r, c) * means(r, c);
            res(r, c) /= (group_sizes[c] - 1);
        }
    }
    
    res = log(res / means);
    res.replace(datum::nan, 0);
    
    return(res);
}


