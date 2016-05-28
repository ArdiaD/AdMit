/* Function which computes the value of the objective function 
   and the gradient of the objective function for the optimization of the
   probabilities 
*/

#include <math.h>

void fnlnf_C(double *lnp, double *lnk, double *lnd, int *Np, int *H, double *f, double *grad)
{
  double tmp, tmp1, tmp2, tmpf1, tmpf2;
  int m, n;
  double f1, f2;
  double lnD[(*Np)*(*H)*(*H)];
  double tmplnd[(*Np)*(*H)], lnw[(*Np)*(*H)];
  double tmpdf1[*H], tmpdf2[*H], df1[*H], df2[*H];

  /* initialization for the objective function */
  m = 0;
  n = 0;
  tmpf1 = 0;
  tmpf2 = 0;
  for (int i=0;i<*Np;i++)
    {
      for (int j=0;j<*H;j++)
	{
	  tmp = 0;
	  for (int k=0;k<*H;k++)
	    {
	      tmp = tmp + exp(lnp[k]+lnd[m]);
	      lnD[m] = lnd[i*(*H)*(*H)+j+k*(*H)];
	      m = m + 1;
	    }
	  tmplnd[n] = log(tmp);
	  lnw[n] = lnk[n]-log(tmp);
	  tmpf1 = tmpf1 + exp(lnp[j]+2*lnw[n]);
	  tmpf2 = tmpf2 + exp(lnp[j]+lnw[n]);  	 	  
	  n = n + 1;
	}
    }
   
  /* initilization for gradient */
  m = 0;  
  n = 0;
  for (int j=0;j<*H;j++)
    {
      tmpdf1[j] = 0;
      tmpdf2[j] = 0;
      df1[j] = 0;
      df2[j] = 0;
    }

  for (int i=0;i<*Np;i++)
    {
      for (int j=0;j<*H;j++)
	{
	  tmp1 = 0;
	  tmp2 = 0;
	  for (int k=0;k<*H;k++)
	    {
	      m = i*(*H)+k;
	      tmp1 = tmp1 + exp(lnp[k]+2*lnk[m]+lnD[n]-3*tmplnd[m]);
	      tmp2 = tmp2 + exp(lnp[k]+lnk[m]+lnD[n]-2*tmplnd[m]);
	      n = n + 1;
	    }
	  tmpdf1[j] = tmpdf1[j] + tmp1;
	  tmpdf2[j] = tmpdf2[j] + tmp2;
	  df1[j] = df1[j] + exp(2*lnw[i*(*H)+j]);
	  df2[j] = df2[j] + exp(lnw[i*(*H)+j]);
	}
    }

  f1 = tmpf1 / (*Np);
  f2 = tmpf2 / (*Np);
  *f = log(f1)-2*log(f2); // output objective function
  
  for (int j=0;j<*H;j++)
    {
      df1[j] = (df1[j] - 2*tmpdf1[j]) / *Np;
      df2[j] = (df2[j] - tmpdf2[j]) / *Np;      
      grad[j] = (df1[j] / f1) - 2 * (df2[j] / f2); //output gradient function
    }
}
