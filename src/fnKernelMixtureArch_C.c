/* Function which computes the kernel of the Mixture of ARCH(1) model */

#include <math.h>

void fnKernelMixtureArch_C(double *theta, int*N, double *y, int *n, double *prior, double *d)
{
  double alpha01, alpha02, alpha1, p, sigma1, sigma2, tmp1, tmp2, tmp;
  int m1, m2;
  for (int j=0;j<*N;j++)
    { /* iterate over parameters */
      m1 = 4*j;
      m2 = 2*j;
      alpha01 = theta[m1]; 
      alpha02 = theta[m1+1];
      alpha1 = theta[m1+2];
      p = theta[m1+3];
      if (prior[m2]==1)
	{ /* if prior satisfied */
	  tmp = 0;
	  for (int i=1;i<*n;i++)
	    { /* compute log-likelihood */
	      sigma1 = sqrt(alpha01 + alpha1 * pow(y[i-1],2));
	      sigma2 = sqrt(alpha02 + alpha1 * pow(y[i-1],2));
	      tmp1 = -0.5 * pow(y[i]/sigma1,2) - log(sigma1);
	      tmp2 = -0.5 * pow(y[i]/sigma2,2) - log(sigma2);
	      tmp = tmp + log( exp(log(p)+tmp1) + exp(log(1-p)+tmp2) );
	    }
	  d[j] = tmp + prior[m2+1]; /* compute kernel */
	}
      else
	{ /* if prior not satisfied */
	  d[j] = prior[m2+1];
	}
    }
}

