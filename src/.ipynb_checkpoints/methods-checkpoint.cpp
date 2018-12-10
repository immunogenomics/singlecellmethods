#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//typedef arma::mat MATTYPE;
//typedef arma::vec VECTYPE;
//typedef arma::fmat MATTYPE;
//typedef arma::fvec VECTYPE;



// [[Rcpp::export]]
arma::mat exp_mean(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, const arma::uvec& groups, const arma::uvec& group_sizes) {
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
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
arma::mat log_vmr(const arma::vec& x, const arma::vec& p, const arma::vec& i, 
                  int ncol, int nrow, const arma::mat& means,
                  const arma::uvec& groups, const arma::uvec& group_sizes) {
    
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
    arma::mat nnzero = arma::zeros<arma::mat>(nrow, ngroups);
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
    res.replace(arma::datum::nan, 0);
    
    return(res);
}

// [[Rcpp::export]]
arma::vec normalizeCLR_dgc(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, int margin) {    
    arma::vec res = x;
    if (margin == 1) {
        // first compute scaling factors for each row
        arma::vec geo_mean = arma::zeros<arma::vec>(nrow);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                // i[j] gives the row num
                geo_mean(i[j]) += std::log1p(x[j]);
            }
        }
        for (int i = 0; i < nrow; i++) {
//            geo_mean(i) = (geo_mean(i) / (1 + ncol));    
            geo_mean(i) = std::exp(geo_mean(i) / ncol);    
        }
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(i[j]));
            }
        }        
    } else {
        // first compute scaling factors for each column
        arma::vec geo_mean = arma::zeros<arma::vec>(ncol);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                geo_mean(c) += std::log1p(x[j]);
            }
            geo_mean(c) = std::exp(geo_mean(c) / nrow);
        }
        
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(c));
            }
        }        
        
    }
    
    return res;
}

/*
// [[Rcpp::export]]
float Hbeta(arma::mat& D, float beta, arma::vec& P, int idx) {
  P = arma::exp(-D.col(idx) * beta);
  float sumP = sum(P);
  float H;
  if (sumP == 0){
      H = 0;
      P = D.col(idx) * 0;
  } else {
      H = log(sumP) + beta * sum(D.col(idx) % P) / sumP;
      P /= sumP;
  }
  return(H);
}

// [[Rcpp::export]]
arma::vec compute_simpson_index(arma::mat& D, arma::umat& knn_idx, arma::vec& batch_labels, int n_batches,
                float perplexity = 15, float tol = 1e-5) {
  int n = D.n_cols;
  arma::vec P = zeros<arma::vec>(D.n_rows);
  arma::vec simpson = zeros<arma::vec>(n);
  float logU = log(perplexity);

  float hbeta, beta, betamin, betamax, H, Hdiff;
  int tries;
  for (int i = 0; i < n ; i++) {
    beta = 1;
    betamin = -datum::inf;
    betamax = datum::inf;
    H = Hbeta(D, beta, P, i);
    Hdiff = H - logU;
    tries = 0;
    // first get neighbor probabilities
    while(std::abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0){
        betamin = beta;
        if (!is_finite(betamax)) beta *= 2;
        else beta = (beta + betamax) / 2;
      } else{
        betamax = beta;
        if (!is_finite(betamin)) beta /= 2;
        else beta = (beta + betamin) / 2;
      }

      H = Hbeta(D, beta, P, i);
      Hdiff = H - logU;
      tries++;
    }

    if (H == 0) {
      simpson.row(i) = -1;
      continue;
    }
    
    // then compute Simpson's Index
    for (int b = 0; b < n_batches; b++) {
      uvec q = find(batch_labels.elem(knn_idx.col(i)) == b); // indices of cells belonging to batch (b)
      if (q.n_elem > 0) {
        float sumP = sum(P.elem(q));
        simpson.row(i) += sumP * sumP;         
      }
    }
  }
  return(simpson);
}








*/