# include <stdio.h>
# include "utils.h"

double *compute_fracs(double *x, int nx){
  int i;
  double *res;
  res = (double *)malloc(sizeof(double) * nalpha);
  for (i=0; i<nx; i++){
	res[i] = 
  }
  
}

double *compute_fracs(NumericMatrix x){
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
	row = dirichlet(REAL(x(i,_)), n);
	fracs(i,_) = row;
  }
  return fracs;
}
