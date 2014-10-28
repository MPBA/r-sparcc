#include <math.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rcpp.h>
#include "utils.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace Rcpp;

/* Extract from dirichlet distribution using gsl library */
double *mydirichlet (double *alpha, int nalpha, int seed){
  int i;
  double rg, fill;
  double *res;
  gsl_rng * r;
  const gsl_rng_type * T;

  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  // Set the random seed
  gsl_rng_set (r, seed);
  
  // Allocate space for result vector
  res = (double *)malloc(sizeof(double) * nalpha);
  
  // fill the result vector
  fill = 1.0/nalpha;
  for (i=0; i<nalpha; i++){
	res[i] = fill;
  }
  
  // extract from gamma distribution
  fill = 0.0;
  for (i=0; i<nalpha; i++){
	// NB add + 1 to input as described in the paper
	rg = gsl_ran_gamma(r, alpha[i] + 1.0, 1.0);
	fill += rg;
	res[i] = rg;
  }
  
  // Normalize
  for (i=0; i<nalpha; i++){
	res[i] /= fill;
  }

  // Free random generator
  gsl_rng_free(r);
  
  return res;
}


/* Function to compute mean of a vector of double */
double mean(double *x, int nx){
  int i;
  double sum=0.0;
  
  // Compute mean
  for (i=0; i<nx; i++){
	sum += x[i];
  }
  return sum/nx;
}


/* Function which compute variance of a double vector */
double var(double *x, int nx){
  int i;
  double meanx, v_tmp;
  
  // Compute mean
  meanx = mean(x, nx);
  
  // Compute variance
  v_tmp = 0.0;
  for (i=0; i<nx; i++){
	v_tmp += pow((x[i] - meanx),2);
  }
  
  // Return
  return v_tmp / (nx- 1);
}

// [[Rcpp::export]]

NumericVector cpp_var (NumericMatrix x){

  int i, j, ncol, nrow;
  double *x_col;
  
  ncol = x.ncol();
  nrow = x.nrow();
  
  NumericVector V(ncol);
  
  x_col = (double *) Calloc(nrow, double); 
  for (i=0; i<ncol; i++){
	
	//fill a double vector from the matrix
	for (j=0; j<nrow; j++){
	  x_col[j] = x(j,i);
	}
	V(i) = var(x_col, nrow);
  }
  
  Free(x_col);
  return V;
}


/* Subtract 2 vectors elemen by element */
double *sub_vect (double *x, double *y, int nx){
  int i;
  double *res;
  
  res = (double *) malloc(nx * sizeof(double));
  
  for (i=0; i<nx; i++){
	res[i] = x[i] - y[i];
  }
  
  return res;
}

/* Compute the variation matrix given T_i,j*/
double  **variation_mat(double **x, int nrow, int ncol){
  
  int i, j;
  double tmp;
  double *mycol;
  double **V;
  
  
  // Alloc space for variation matrix (p x p)
  V = (double **) malloc(ncol * sizeof(double *));
  for (i=0; i<ncol; i++){
	V[i] = (double *)malloc(ncol * sizeof(double));
  }
  
  V[0][0] = 1.0;
  //  Compute variation matrix
  for (i=0; i<(ncol-1); i++){
	for (j=i+1; j<ncol; j++){
	  mycol = sub_vect(x[i], x[j], nrow);
	  tmp = var(mycol, nrow);
	  V[i][j] = tmp;
	  V[j][i] = tmp;
	  V[j][j] = 1.0;
	}
  }
  
  return V;
}

/* Convert conts to relative fractions */
double *counts2frac(double *x, int nx, int seed){
  int j;
  double *res, *tmp;
  
  res = (double *) malloc(nx * sizeof(double));
  
  // Extract from a dirichlet distribution
  tmp = mydirichlet(x, nx, seed);
  
  // Put the results into the FRACS matrix
  for(j=0; j<nx; j++){
	 res[j]= log(tmp[j]);
  }
  
  free(tmp);
  return res;
}

/* invert matrix */
gsl_matrix *inv_mat(double **x, int nrow, int ncol){
  int i, j, s;
  gsl_matrix *M = gsl_matrix_alloc(nrow, ncol);
  gsl_matrix *Minv = gsl_matrix_alloc(nrow, ncol);

  for (i=0; i<ncol; i++){
	for(j=0; j<nrow; j++){
	  gsl_matrix_set(M, j, i, x[i][j]);
	}
  }
  gsl_permutation * p = gsl_permutation_alloc (nrow);

  gsl_linalg_LU_decomp (M, p, &s);
  gsl_linalg_LU_invert (M, p, Minv);

  gsl_permutation_free (p);
  gsl_matrix_free (M);
  
  return Minv;
}

/* sum over columns */
double *colsum (double **x, int nrow, int ncol){
  int i, j;
  double *sum;
  
  sum = (double *)malloc(ncol * sizeof(double));
  
  for (i=0; i<ncol; i++){
	sum[i] = 0.0;
	  for (j=0; j<nrow; j++){
		sum[i] += x[i][j];
	  }
  }
  return sum;
}

/* Put here the function that estimates the basis variances */
/* Need to invert a matrix */
double estimate_basis_variance(double **V, int nx, double V_min, int *excl_idx, gsl_matrix *Cov_mat, gsl_matrix *M){
  
  
  
}




/* Convert a R matrix to a double pointer' array in C */
double **rmat2double(NumericMatrix x, int nrow, int ncol){
  int i, j;
  double **res;
  
  // Allocate space for the copy of the matrix
  res = (double **) Calloc(ncol, double *);
  for (i=0; i<ncol; i++){
	res[i] = (double *) Calloc(nrow, double);
	for (j=0; j<nrow; j++){
	  res[i][j] = x(j,i);
	}
  }
  
  return res;
}

// [[Rcpp::export]]

NumericMatrix sparcc_all(NumericMatrix x){
  int i, j, nrow, ncol;
  double *myrow, *tmp;
  double **counts, **fracs, **V;
  
  nrow = x.nrow();
  ncol = x.ncol();
  
  NumericMatrix Vmat(ncol, ncol);

  // Alloc matrix relative fractions
  fracs = (double **) Calloc (ncol, double *);
  for (i=0; i<ncol; i++){
  	fracs[i] = (double *) Calloc (nrow, double);
  }
  
  // Alloc row vector
  // row = (double)malloc(ncol*sizeof(double));
  myrow = (double *) Calloc(ncol, double);
  
  //Copy R matrix
  counts = rmat2double(x, nrow, ncol);
  
  for (i=0; i<nrow; i++){
  	for(j=0; j<ncol; j++){
  	  myrow[j] = counts[j][i];
  	}
	tmp = counts2frac(myrow, ncol, i);
	for(j=0; j<ncol; j++){
	  fracs[j][i] = tmp[j];
	}
  }
  
  V = variation_mat(fracs, nrow, ncol);
  
  Vmat(0,0) = V[0][0];
  // Copy variation matrix results to Vmat
  for (i=0; i<ncol-1; i++){
  	for (j=i+1; j<ncol; j++){
  	  Vmat(i,j) = V[i][j];
  	  Vmat(j,i) = V[j][i];
  	  Vmat(j,j) = V[j][j];
  	}
  }
  
  //free counts
  for (i=0; i<ncol; i++){
  	free(counts[i]);
  }
  free(counts);
  
  //free fracs
  for (i=0; i<ncol; i++){
	Free(fracs[i]);
  }
  Free(fracs);
  
  // free V
  for (i=0; i<ncol; i++){
  	free(V[i]);
  }
  free(V);
  
  Free(myrow);
  
  return Vmat;
}


