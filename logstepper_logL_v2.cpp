#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

const double log2pi2 = log(2.0 * M_PI)/2.0;

// [[Rcpp::export]]
double logstep(arma::mat X, arma::vec weights, arma::mat Mean, arma::mat Sigma,
               arma::vec & A0, arma::mat & A1, arma::mat & A2) {
  
  int n = X.n_rows;
  int p = Sigma.n_cols;
  int K = Mean.n_rows;
  double constant = -p*log2pi2;
  arma::mat Rinv(K*p,p); 
  arma::vec logSqrtDetvarcovM(K);
  for (int k = 0; k < K; ++k) {
    arma::mat S = Sigma.submat(k*p,0,k*p + (p-1),p-1);
    Rinv.submat(k*p,0,k*p + (p-1),p-1) = inv(trimatu(chol(S)));
    logSqrtDetvarcovM[k] = sum(log(Rinv.diag()));
  }
  double b;
  double a;
  double logL = 0;
  
  for (int i = 0; i < n; ++i) {
    
    arma::colvec y = trans(X.row(i));
    arma::vec logprobs(K);
    arma::vec probs(K);
    
    a = 0;
    b = 0;
    
    for (int k = 0; k < K; ++k) {
      
      arma::colvec M = trans(Mean.row(k));
      
      arma::colvec x_i = y - M;
      arma::rowvec xRinv = trans(x_i) * Rinv.submat(k*p,0,k*p + (p-1),p-1);
      
      double quadform = sum(xRinv%xRinv);
      //double GAUSS = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
      double lognorm = -0.5*quadform + logSqrtDetvarcovM[k] + constant;
      
      //probs[k] = weights[k] * GAUSS;
      if (k<1) {
        a = log(weights[k]) + lognorm;
        logprobs[k] = a;
      } else {
        b = log(weights[k]) + lognorm;
        if (a >= b) {
          a += log(1 + exp(b - a));
        } else {
          a = b + log(1 + exp(a - b));
        }
        logprobs[k] = b;
      }
    }
    //double SUM = sum(probs);
    
    logL += a;
    
    //probs /= SUM;
    for (int k = 0; k < K; ++k) {
      probs[k] = exp(logprobs[k] - a);
      
      A1.row(k) += probs[k] * arma::trans(y);
      A2.submat(k*p,0,k*p + (p-1),p-1) += probs[k] * (y * arma::trans(y));
    }
    A0 += probs;
  }
  
  return logL;
}

