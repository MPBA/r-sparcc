#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

int main (int argc, char *argv[])
{
  /* set up GSL RNG */
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  /* end of GSL setup */

  int i,n;
  double gauss,gamma;  

  n=atoi(argv[1]);
  for (i=0;i<n;i++)
    {
      gauss=gsl_ran_gaussian(r,2.0);
      gamma=gsl_ran_gamma(r,2.0,3.0);
      printf("%2.4f %2.4f\n", gauss,gamma);
    }
  return(0);
}
