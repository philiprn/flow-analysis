#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

const double log2pi2 = log(2.0 * M_PI)/2.0;

// [[Rcpp::export()]]
arma::ucolvec getClusters(arma::mat X, arma::vec weights, arma::mat Mean, arma::mat Sigma) {
  
  int n = X.n_rows;
  int p = Sigma.n_cols;
  int K = Mean.n_rows;
  //arma::vec clusterlabels(n);
  
  double constant = -p*log2pi2;
  
  arma::mat logprobs(n,K);
  
  for (int k = 0; k < K; ++k) {
    arma::mat S = Sigma.submat(k*p,0,k*p + (p-1),p-1);
    arma::mat Rinv = inv(trimatu(chol(S)));
    double logSqrtDetvarcovM = sum(log(Rinv.diag()));
    arma::colvec M = trans(Mean.row(k));
    for (int i = 0; i < n; ++i) {
      
      arma::colvec y = trans(X.row(i));
      
      arma::colvec x_i = y - M;
      arma::rowvec xRinv = trans(x_i)*Rinv;
      
      double quadform = sum(xRinv%xRinv);
      double lognorm = -0.5*quadform + logSqrtDetvarcovM + constant;
      
      logprobs(i,k) = lognorm + log(weights[k]);
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  
  arma::ucolvec clusterlabels = index_max(logprobs,1) + 1;
  
  return clusterlabels;
}
