#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::rowvec getm(arma::mat X) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  arma::rowvec means(p);
  
  for (int i = 0; i < n; ++i) {
    means += X.row(i);
  }
  
  means /= n;
  
  return means;
}

