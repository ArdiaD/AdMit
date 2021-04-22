/* Function which generates the MH chain */

#include <math.h>

void fnMH_C(double *theta, int *N, int *k, double *lnw, double *u, double *draws, int *ns)
{
  int m, s;
  double e;
  double tmp[*k];

  /* initilize chain */
  for (int j=0;j<*k;j++)
    { 
      tmp[j] = theta[j];
      draws[j] = theta[j];
    }
  
  s = 0;
  /* iterate over the chain */
  for (int i=1;i<*N;i++)
    {
      e = exp(lnw[i]-lnw[s]);
      if (e>1)
	{ /* minimum between e and 1 */
	  e = 1;
	}
     
      m = i*(*k);
      /* compare with random U(0,1) */
      if (u[i] <= e)
	{ 
	  for (int j=0;j<*k;j++)
	    {
	      tmp[j] = theta[m+j];
	    }
	  s = i;           /* new past position */
	  *ns = *ns+1;     /* increment counter for acceptance rate */
	}
      /* assign tmp to draws */
      for (int j=0;j<*k;j++)
	draws[m+j] = tmp[j];
    }
}
