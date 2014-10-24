#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericVector myrdirichlet(NumericVector alpha){
  int l = alpha.size();
  double sm;
  NumericVector x(l), res(l);
  NumericVector::iterator it;
  for (int i=0; i<l; i++){
	x(i) = rgamma(1, alpha(i))[0];
  }
  sm = std::accumulate(x.begin(), x.end(), 0.0);
  for (it = x.begin(); it != x.end(); ++it){
	*it /= sm;
  }
  return x;
}

// [[Rcpp::export]]

NumericMatrix compute_fracs(NumericMatrix x){
  /* Convert count matrix to fractions */
  int n=x.nrow();
  int p=x.ncol();
  double initval = 1/p;
  NumericVector row, tmp;
  
  NumericMatrix fracs(n,p);
  NumericMatrix::iterator it;
  
  /* Initialize fracs matrix*/
  for (it = fracs.begin(); it != fracs.end(); ++it){
	*it = initval;
  }
  
  /* Sum 1 to x */
  for (it = x.begin(); it != x.end(); ++it){
  	*it += 1.0;
  }
  
  /* Sample from dirichlet */
  for (int i=0; i<n; i++){
	row = myrdirichlet(x(i,_));
	fracs(i,_) = row;
  }
  return fracs;
}
