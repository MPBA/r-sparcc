#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double *mydirichlet (double *alpha, int nalpha){
  int i;
  double rg;
  double *res;
  gsl_rng * r;
  const gsl_rng_type * T;
  
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  res = (double *)malloc(sizeof(double) * nalpha);
  
  for (i=0; i<nalpha; i++){
	// NB add + 1 to input as described in the paper
	rg = gsl_ran_gamma(r, alpha[i] + 1, 1.0);
	res[i] = rg;
  }
  
  gsl_rng_free(r);
  
  return res;
}




int main(){
  double *alpha, *res;
  int nalpha = 4;
  int i;
  
  alpha = (double *)malloc(sizeof(double) * nalpha);
  
  for (i = 0; i<nalpha; i++){
  	alpha[i] = i*5.0;
  }
  res = mydirichlet(alpha, nalpha);
  
  for (i = 0; i<nalpha; i++){
  	printf("res[%d]=%.3f \n", i, res[i]);
  	//printf(i);
  }
  
  //free results
  free(res);
  
  return 0;
}
