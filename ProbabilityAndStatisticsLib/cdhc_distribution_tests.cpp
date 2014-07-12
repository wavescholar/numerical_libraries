 /*******************************
 * Copyright (c) <2007>, <Bruce Campbell> All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  *  
 * Bruce B Campbell 07 08 2014  *
 ********************************/
//Distribution tests.
//GRASS GIS distributes this code as well, I have a note somewhere with the true origin.    

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int dcmp(const void* ,const void*);
double normp (double z);
double alnfac (int j);
double alnorm (double x,int  upper);
double correc (int i,int  n);
double ppnd16 (double p);
double poly (double c[], int nord, double x);
double xinormal (double pee);
double * coeff_variation(double * x,int n);
double *mod_maxlik_ratio (double * x,int  n);

int dcmp (const void* i,const void*  j)
   
{
	void* b=const_cast< void*>(i);
	double* l=reinterpret_cast< double*>(&b);
//	double rt=*r;
	void* c=const_cast< void*>(j);
	double* r=reinterpret_cast< double*>(&c);

  if (*l<*r)
        return  -1;
 
    if (*l >*r)
        return 1;
 
    return 0;
}

double *anderson_darling_exp (double* x,int  n)
  {
  static double y[2];
  double sqrt2, mean = 0.0, *xcopy, fx, sum3 = 0.0;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in anderson_darling\n"), exit (-1);

  sqrt2 = sqrt ((double) 2.0);

  for (i = 0; i < n; ++i)
  {
    xcopy[i]= x[i];
    mean += x[i];
  }
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    fx = 1-exp(-xcopy[i]/mean);
    sum3 += (2.0 * i + 1) * (log (fx) -xcopy[n-i-1]/mean);
  }
  y[0] = (1.0+0.3/n)*(-n-sum3/n);
#ifdef NOISY
  printf ("  TEST20 AD(E)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  free(xcopy);
  return y;
}

double *anderson_darling (double* x,int  n)
  {
  int i;
  static double y[2];
  double sqrt2, mean = 0.0, sdx=0.0,*xcopy, fx;
  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in anderson_darling\n"), exit (-1);

  sqrt2 = sqrt ((double) 2.0);
  y[0]=y[1]=0.0;

  for (i = 0; i < n; ++i)
  {
    xcopy[i]= x[i];
    mean += x[i];
    sdx += x[i] * x[i];
  }
  sdx = sqrt ((n * sdx - mean * mean) / (n * (n - 1.0)));
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
    xcopy[i] = (xcopy[i] - mean) / sdx;
  for (i = 0; i < n; ++i)
  {
    fx = 0.5 + normp (xcopy[i] / sqrt2) / 2.0;
    if (fx <= 1e-5)
      fx = 1e-5;
    if (fx >= .99999)
      fx = 0.99999;
    y[1] += (2.0*i+1.0)  * log(fx) + (2.0*(n-i)-1.0) * log(1-fx);
  }
  y[1]=-n-y[1]/n;
  y[0] = y[1] * (0.75 / n + 1.0 + 2.25 / (n * n));

#ifdef NOISY
  printf ("  TEST8  AD(N)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  free(xcopy);
  return y;
}
/*-Algorithm AS 177
 * Expected Normal Order Statistics (Exact and Approximate),
 * by J.P. Royston, 1982.
 * Applied Statistics, 31(2):161-165.
 *
 * Translation to C by James Darrell McCauley, mccauley@ecn.purdue.edu.
 *
 * The functions nscor1() and nscor2() calculate the expected values of
 * normal order statistics in exact or approximate form, respectively.
 *
 */

#define NSTEP 721
#define H 0.025

void nscor1 (double s[], int n, int n2, double work[], int* ifault)

/* exact calculation of normal scores */
{
  double ani, c, c1, d, scor;
  int i,j;

  *ifault = 3;
  if (n2 != n / 2)
    return;
  *ifault = 1;
  if (n <= 1)
    return;
  *ifault = 0;
  if (n > 2000)
    *ifault = 2;
  /* calculate the natural log of factorial(n) */
  c1 = alnfac (n);
  d = c1 - log ((double) n);
  /* accumulate ordinates for calculation of integral for rankits */
  for (i = 0; i < n2; ++i)
  {
    ani = (double) n - i - 1;
    c = c1 - d;
    for (scor=0.0, j = 0; j < NSTEP; ++j)
      scor += work[0*NSTEP+j] *
	exp (work[1*NSTEP+j] + work[2*NSTEP+j] * i
	     + work[3*NSTEP+j] * ani + c);
    s[i] = scor * H;
    d += log ((double)(i + 1.0) / ani);
  }
  return;
}

void init (double work[])
  
{
  double xstart = -9.0, pi2 = -0.918938533, xx;
  int i;

  xx = xstart;
  /* set up arrays for calculation of integral */
  for (i = 0; i < NSTEP; ++i)
  {
    work[0*NSTEP+i] = xx;
    work[1*NSTEP+i] = pi2 - xx * xx * 0.5;
    work[2*NSTEP+i] = log (alnorm (xx,1));
    work[3*NSTEP+i] = log (alnorm (xx,0));
    xx = xstart + H * (i + 1.0);
  }
  return;
}

double alnfac (int j)
 
/*-Algorithm AS 177.2 Appl. Statist. (1982) Vol.31, No.2
 * Natural logarithm of factorial for non-negative argument
 */
{
  static double r[7] = {0.0, 0.0, 0.69314718056, 1.79175946923,
  3.17805383035, 4.78749174278, 6.57925121101};
  double w, z;

  if (j == 1)
    return (double) 1.0;
  else if (j <= 7)
    return r[j];
  w = (double) j + 1;
  z = 1.0 / (w * w);
    return (w - 0.5) * log (w) - w + 0.918938522305 +
    (((4.0 - 3.0 * z) * z - 14.0) * z + 420.0) / (5040.0 * w);
}

void nscor2 (double s[], int n, int n2,int * ifault)
 
/*-Algorithm AS 177.3 Appl. Statist. (1982) Vol.31, No.2
 * Approximation for Rankits
 */
{
  static double eps[4] = {0.419885, 0.450536, 0.456936, 0.468488};
  static double dl1[4] = {0.112063, 0.121770, 0.239299, 0.215159};
  static double dl2[4] = {0.080122, 0.111348, -0.211867, -0.115049};
  static double gam[4] = {0.474798, 0.469051, 0.208597, 0.259784};
  static double lam[4] = {0.282765, 0.304856, 0.407708, 0.414093};
  static double bb = -0.283833, d = -0.106136, b1 = 0.5641896;
  double e1, e2, l1;
  int i, k;

  *ifault = 3;
  if (n2 != n / 2)
    return;
  *ifault = 1;
  if (n <= 1)
    return;
  *ifault = 0;
  if (n > 2000)
    *ifault = 2;
  s[0] = b1;
  if (n == 2)
    return;
  /* calculate normal areas for 3 largest rankits */
  k=(n2<3) ? n2 : 3;
  for (i = 0; i < k; ++i)
  {
    e1 = (1.0 + i - eps[i]) / (n + gam[i]);
    e2 = pow (e1, lam[i]);
    s[i] = e1 + e2 * (dl1[i] + e2 * dl2[i]) / n - correc (1 + i, n);
  }
  if (n2 != k)
  {
    /* calculate normal areas for remaining rankits */
    for (i = 3; i < n2; ++i)
    {
      l1 = lam[3] + bb / (1.0 + i + d);
      e1 = (1.0 + i - eps[3]) / (n + gam[3]);
      e2 = pow (e1, l1);
      s[i] = e1 + e2 * (dl1[3] + e2 * dl2[3]) / n - correc (1+i, n);
    }
  }
  /* convert normal tail areas to normal deviates */
  for (i = 0; i < n2; ++i)
    s[i] = -ppnd16 (s[i]);
  return;
}

double correc (int i,int  n)
 /*-Algorithm AS 177.4 Appl. Statist. (1982) Vol.31, No.2
 * Calculates correction for tail area of noraml distribution
 * corresponding to ith largest rankit in sample size n.
 */
{
  static double c1[7] = {9.5, 28.7, 1.9, 0.0, -7.0, -6.2, -1.6};
  static double c2[7] = {-6.195e3, -9.569e3, -6.728e3, -17.614e3,
  -8.278e3, -3.570e3, 1.075e3};
  static double c3[7] = {9.338e4, 1.7516e5, 4.1040e5, 2.157e6,
  2.376e6, 2.065e6, 2.065e6};
  static double mic = 1.0e-6, c14 = 1.9e-5;
  double an, ret_val;

  ret_val = c14;
  if (i * n == 4)
    return ret_val;
  ret_val = 0.0;
  if (i < 1 || i > 7)
    return ret_val;
  else if (i != 4 && n > 20)
    return ret_val;
  else if (i == 4 && n > 40)
    return ret_val;
  else
  {
    an = 1.0 / (double) (n * n);
    ret_val = (c1[i-1] + an * (c2[i-1] + an * c3[i-1])) * mic;
    return ret_val;
  }
}
//
//void wext_as181 (x, n, ssq, a, n2, eps, w, pw, ifault)
//  double x[], ssq, a[], eps, *w, *pw;
//  int n, n2, *ifault;
//
///*-Algorithm AS 181
// * by J.P. Royston, 1982.
// * Applied Statistics 31(2):176-180
// *
// * Translation to C by James Darrell McCauley, mccauley@ecn.purdue.edu.
// *
// * Calculates Shapiro and Wilk's W statistic and its sig. level
// *
// * Originally used:
// * Auxiliary routines required: ALNORM = algorithm AS 66 and NSCOR2
// * from AS 177.
//
// * Note: ppnd() from as66 was replaced with ppnd16() from as241.
// */
//{
//  double eu3, lamda, ybar, sdy, al, un, ww, y, z;
//  int i, j, n3, nc;
//  static double wa[3] = {0.118898, 0.133414, 0.327907};
//  static double wb[4] = {-0.37542, -0.492145, -1.124332, -0.199422};
//  static double wc[4] = {-3.15805, 0.729399, 3.01855, 1.558776};
//  static double wd[6] = {0.480385, 0.318828, 0.0, -0.0241665, 0.00879701,
//  0.002989646};
//  static double we[6] = {-1.91487, -1.37888, -0.04183209, 0.1066339,
//  -0.03513666, -0.01504614};
//  static double wf[7] = {-3.73538, -1.015807, -0.331885, 0.1773538,
//  -0.01638782, -0.03215018, 0.003852646};
//  static double unl[3] = {-3.8, -3.0, -1.0};
//  static double unh[3] = {8.6, 5.8, 5.4};
//  static int nc1[3] = {5, 5, 5};
//  static int nc2[3] = {3, 4, 5};
//  double c[5];
//  int upper = 1;
//  static double pi6 = 1.90985932, stqr = 1.04719755;
//  static double zero = 0.0, tqr = 0.75, one = 1.0;
//  static double onept4 = 1.4, three = 3.0, five = 5.0;
//  static double c1[5][3] = {
//    {-1.26233, -2.28135, -3.30623},
//    {1.87969, 2.26186, 2.76287},
//    {0.0649583, 0.0, -0.83484},
//    {-0.0475604, 0.0, 1.20857},
//    {-0.0139682, -0.00865763, -0.507590}
//  };
//  static double c2[5][3] = {
//    {-0.287696, -1.63638, -5.991908},
//    {1.78953, 5.60924, 21.04575},
//    {-0.180114, -3.63738, -24.58061},
//    {0.0, 1.08439, 13.78661},
//    {0.0, 0.0, -2.835295}
//  };
//  double poly(), alnorm();
//
//  *ifault = 1;
//
//  *pw = one;
//  *w = one;
//  if (n <= 2)
//    return;
//  *ifault = 3;
//  if (n / 2 != n2)
//    return;
//  *ifault = 2;
//  if (n > 2000)
//    return;
//  *ifault = 0;
//  i = n - 1;
//  for (*w = 0.0, j = 0; j < n2; ++j)
//    *w += a[j] * (x[i--] - x[j]);
//  *w *= *w / ssq;
//  if (*w > one)
//  {
//    *w = one;
//    return;
//  }
//  else if (n > 6)		/* Get significance level of W */
//  {
//    /*
//     * N between 7 and 2000 ... Transform W to Y, get mean and sd,
//     * standardize and get significance level
//     */
//
//    if (n <= 20)
//    {
//      al = log ((double) n) - three;
//      lamda = poly (wa, 3, al);
//      ybar = exp (poly (wb, 4, al));
//      sdy = exp (poly (wc, 4, al));
//    }
//    else
//    {
//      al = log ((double) n) - five;
//      lamda = poly (wd, 6, al);
//      ybar = exp (poly (we, 6, al));
//      sdy = exp (poly (wf, 7, al));
//    }
//
//    y = pow (one - *w, lamda);
//    z = (y - ybar) / sdy;
//    *pw = alnorm (z, upper);
//    return;
//  }
//  else
//  {
//    /* Deal with N less than 7 (Exact significance level for N = 3). */
//    if (*w >= eps)
//    {
//      ww = *w;
//      if (*w >= eps)
//      {
//	ww = *w;
//	if (n == 3)
//	{
//	  *pw = pi6 * (atan (sqrt (ww / (one - ww))) - stqr);
//	  return;
//	}
//	un = log ((*w - eps) / (one - *w));
//	n3 = n - 3;
//	if (un >= unl[n3 - 1])
//	{
//	  if (un <= onept4)
//	  {
//	    nc = nc1[n3 - 1];
//	    for (i = 0; i < nc; ++i)
//	      c[i] = c1[i][n3 - 1];
//	    eu3 = exp (poly (c, nc, un));
//	  }
//	  else
//	  {
//	    if (un > unh[n3 - 1])
//	      return;
//	    nc = nc2[n3 - 1];
//	    for (i = 0; i < nc; ++i)
//	      c[i] = c2[i][n3 - 1];
//
//	    un = log (un); /*alog*/
//	    eu3 = exp (exp (poly (c, nc, un)));
//	  }
//	  ww = (eu3 + tqr) / (one + eu3);
//	  *pw = pi6 * (atan (sqrt (ww / (one - ww))) - stqr);
//	  return;
//	}
//      }
//    }
//    *pw = zero;
//    return;
//  }
//}
void wext (double x[],int n,double ssq, double a[], int n2,double  eps,double*  w,double*  pw,int* ifault)
  
/*-Algorithm AS 181
 * by J.P. Royston, 1982.
 * Applied Statistics 31(2):176-180
 *
 * Translation to C by James Darrell McCauley, mccauley@ecn.purdue.edu.
 *
 * Calculates Shapiro and Wilk's W statistic and its sig. level
 *
 * Originally used:
 * Auxiliary routines required: ALNORM = algorithm AS 66 and NSCOR2
 * from AS 177.

 * Note: ppnd() from as66 was replaced with ppnd16() from as241.
 */
{
  double eu3, lamda, ybar, sdy, al, un, ww, y, z;
  int i, j, n3, nc;
  static double wa[3] = {0.118898, 0.133414, 0.327907};
  static double wb[4] = {-0.37542, -0.492145, -1.124332, -0.199422};
  static double wc[4] = {-3.15805, 0.729399, 3.01855, 1.558776};
  static double wd[6] = {0.480385, 0.318828, 0.0, -0.0241665, 0.00879701,
  0.002989646};
  static double we[6] = {-1.91487, -1.37888, -0.04183209, 0.1066339,
  -0.03513666, -0.01504614};
  static double wf[7] = {-3.73538, -1.015807, -0.331885, 0.1773538,
  -0.01638782, -0.03215018, 0.003852646};
  static double unl[3] = {-3.8, -3.0, -1.0};
  static double unh[3] = {8.6, 5.8, 5.4};
  static int nc1[3] = {5, 5, 5};
  static int nc2[3] = {3, 4, 5};
  double c[5];
  int upper = 1;
  static double pi6 = 1.90985932, stqr = 1.04719755;
  static double zero = 0.0, tqr = 0.75, one = 1.0;
  static double onept4 = 1.4, three = 3.0, five = 5.0;
  static double c1[5][3] = {
    {-1.26233, -2.28135, -3.30623},
    {1.87969, 2.26186, 2.76287},
    {0.0649583, 0.0, -0.83484},
    {-0.0475604, 0.0, 1.20857},
    {-0.0139682, -0.00865763, -0.507590}
  };
  static double c2[5][3] = {
    {-0.287696, -1.63638, -5.991908},
    {1.78953, 5.60924, 21.04575},
    {-0.180114, -3.63738, -24.58061},
    {0.0, 1.08439, 13.78661},
    {0.0, 0.0, -2.835295}
  };


  *ifault = 1;

  *pw = one;
  *w = one;
  if (n <= 2)
    return;
  *ifault = 3;
  if (n / 2 != n2)
    return;
  *ifault = 2;
  if (n > 2000)
    return;
  *ifault = 0;
  i = n - 1;
  for (*w = 0.0, j = 0; j < n2; ++j)
    *w += a[j] * (x[i--] - x[j]);
  *w *= *w / ssq;
  if (*w > one)
  {
    *w = one;
    return;
  }
  else if (n > 6)		/* Get significance level of W */
  {
    /*
     * N between 7 and 2000 ... Transform W to Y, get mean and sd,
     * standardize and get significance level
     */

    if (n <= 20)
    {
      al = log ((double) n) - three;
      lamda = poly (wa, 3, al);
      ybar = exp (poly (wb, 4, al));
      sdy = exp (poly (wc, 4, al));
    }
    else
    {
      al = log ((double) n) - five;
      lamda = poly (wd, 6, al);
      ybar = exp (poly (we, 6, al));
      sdy = exp (poly (wf, 7, al));
    }

    y = pow (one - *w, lamda);
    z = (y - ybar) / sdy;
    *pw = alnorm (z, upper);
    return;
  }
  else
  {
    /* Deal with N less than 7 (Exact significance level for N = 3). */
    if (*w >= eps)
    {
      ww = *w;
      if (*w >= eps)
      {
	ww = *w;
	if (n == 3)
	{
	  *pw = pi6 * (atan (sqrt (ww / (one - ww))) - stqr);
	  return;
	}
	un = log ((*w - eps) / (one - *w));
	n3 = n - 3;
	if (un >= unl[n3 - 1])
	{
	  if (un <= onept4)
	  {
	    nc = nc1[n3 - 1];
	    for (i = 0; i < nc; ++i)
	      c[i] = c1[i][n3 - 1];
	    eu3 = exp (poly (c, nc, un));
	  }
	  else
	  {
	    if (un > unh[n3 - 1])
	      return;
	    nc = nc2[n3 - 1];
	    for (i = 0; i < nc; ++i)
	      c[i] = c2[i][n3 - 1];
	    un = log (un); /*alog*/
	    eu3 = exp (exp (poly (c, nc, un)));
	  }
	  ww = (eu3 + tqr) / (one + eu3);
	  *pw = pi6 * (atan (sqrt (ww / (one - ww))) - stqr);
	  return;
	}
      }
    }
    *pw = zero;
    return;
  }
}

void wcoef (double a[],int  n,int  n2,double* eps,int* ifault)
 
/*
 * Algorithm AS 181.1   Appl. Statist.  (1982) Vol. 31, No. 2
 * 
 * Obtain array A of weights for calculating W
 */
{
  static double c4[2] = {0.6869, 0.1678};
  static double c5[2] = {0.6647, 0.2412};
  static double c6[3] = {0.6431, 0.2806, 0.0875};
  static double rsqrt2 = 0.70710678;
  double a1star, a1sq, sastar, an;
  int j;

  *ifault = 1;
  if (n <= 2)
    return;
  *ifault = 3;
  if (n / 2 != n2)
    return;
  *ifault = 2;
  if (n > 2000)
    return;
  *ifault = 0;
  if (n > 6)
  {
    /* Calculate rankits using approximate function nscor2().  (AS177) */

    nscor2 (a, n, n2, ifault);
    for (sastar = 0.0, j = 1; j < n2; ++j)
      sastar += a[j] * a[j];
    sastar *= 8.0;

    an = n;
    if (n <= 20)
      an--;
    a1sq = exp (log (6.0 * an + 7.0) - log (6.0 * an + 13.0)
		+ 0.5 * (1.0 + (an - 2.0) * log (an + 1.0) - (an - 1.0)
			 * log (an + 2.0)));
    a1star = sastar / (1.0 / a1sq - 2.0);
    sastar = sqrt (sastar + 2.0 * a1star);
    a[0] = sqrt (a1star) / sastar;
    for (j = 1; j < n2; ++j)
      a[j] = 2.0 * a[j] / sastar;
  }
  else
  {
    /* Use exact values for weights */

    a[0] = rsqrt2;
    if (n != 3)
    {
      if (n - 3 == 3)
         for(j=0;j<3;++j) a[j]=c6[j];
      else if (n - 3 == 2)
         for(j=0;j<2;++j) a[j]=c5[j];
      else
         for(j=0;j<2;++j) a[j]=c4[j];
      /*-
            goto (40,50,60), n3
         40 do 45 j = 1,2
         45 a(j) = c4(j)
            goto 70
         50 do 55 j = 1,2
         55 a(j) = c5(j)
            goto 70
         60 do 65 j = 1,3
         65 a(j) = c6(j)
      */
    }
  }

  /* Calculate the minimum possible value of W */
  *eps = a[0] * a[0] / (1.0 - 1.0 / (double) n);
  return;
}

double poly (double c[], int nord, double x)

/*
 * Algorithm AS 181.2   Appl. Statist.  (1982) Vol. 31, No. 2
 * 
 * Calculates the algebraic polynomial of order nored-1 with array of
 * coefficients c.  Zero order coefficient is c(1)
 */
{
  double p;
  int n2, i, j;

  if (nord == 1)
    return c[0];
  p = x * c[nord - 1];
  if (nord != 2)
  {
    n2 = nord - 2;
    j = n2;
    for (i = 0; i < n2; ++i)
      p = (p + c[j--]) * x;
  }
  return c[0] + p;
}

void wgp (double x[],int  n,double ssq, double gp, double h,double  a[],int  n2,double  eps,double  w,double  u,double  p,int*  ifault)
/*
 * AS R63 Appl. Statist. (1986) Vol. 35, No.2
 * 
 * A remark on AS 181
 * 
 * Calculates Sheppard corrected version of W test.
 * 
 * Auxiliary functions required: ALNORM = algorithm AS 66, and PPND =
 * algorithm AS 111 (or PPND7 from AS 241).
 */
{
  double zbar, zsd, an1, hh;

  zbar = 0.0;
  zsd = 1.0;
  *ifault = 1;
  if (n < 7)
    return;
  if (gp > 0.0)			/* No correction applied if gp=0. */
  {
    an1 = (double) (n - 1);
    /* correct ssq and find standardized grouping interval (h) */
    ssq = ssq - an1 * gp * gp / 12.0;
    h = gp / sqrt (ssq / an1);
    *ifault = 4;
    if (h > 1.5)
      return;
  }
  wext (x, n, ssq, a, n2, eps, &w,&p, ifault);
  if (*ifault != 0)
    return;
  if (!(p > 0.0 && p < 1.0))
  {
    u = 5.0 - 10.0 * p;
      return;
  }
  if (gp > 0.0)
  {
    /* correct u for grouping interval (n<=100 and n>100 separately) */
    hh = sqrt (h);
    if (n <= 100)
    {
      zbar = -h * (1.07457 + hh * (-2.8185 + hh * 1.8898));
      zsd = 1.0 + h * (0.50933 + hh * (-0.98305 + hh * 0.7408));
    }
    else
    {
      zbar = -h * (0.96436 + hh * (-2.1300 + hh * 1.3196));
      zsd = 1.0 + h * (0.2579 + h * 0.15225);
    }
  }
  /* ppnd is AS 111 (Beasley and Springer, 1977) */
  u = (-ppnd16 (p) - zbar) / zsd;
  /* alnorm is AS 66 (Hill, 1973) */
  p = alnorm (u,1);
  return;
}

double ppnd7 (double p)

/*-
 * algorithm as241  appl. statist. (1988) 37(3):477-484.
 * produces the normal deviate z corresponding to a given lower tail
 * area of p; z is accurate to about 1 part in 10**7.
 *
 * the hash sums below are the sums of the mantissas of the coefficients.
 * they are included for use in checking transcription.
 */
{
  static double zero = 0.0, one = 1.0, half = 0.5;
  static double split1 = 0.425, split2 = 5.0;
  static double const1 = 0.180625, const2 = 1.6;

  /* coefficients for p close to 0.5 */
  static double a[4] = {3.3871327179, 5.0434271938e+01,
  1.5929113202e+02, 5.9109374720e+01};
  static double b[4] = {0.0, 1.7895169469e+01, 7.8757757664e+01,
  6.7187563600e+01};

  /* hash sum ab    32.3184577772 */
  /* coefficients for p not close to 0, 0.5 or 1. */
  static double c[4] = {1.4234372777e+00, 2.7568153900e+00,
  1.3067284816e+00, 1.7023821103e-01};
  static double d[3] = {0.0, 7.3700164250e-01, 1.2021132975e-01};

  /* hash sum cd    15.7614929821 */
  /* coefficients for p near 0 or 1. */
  static double e[4] = {6.6579051150e+00, 3.0812263860e+00,
  4.2868294337e-01, 1.7337203997e-02};
  static double f[3] = {0.0, 2.4197894225e-01, 1.2258202635e-02};

  /* hash sum ef    19.4052910204 */
  double q, r, ret;

  q = p - half;
  if (fabs (q) <= split1)
  {
    r = const1 - q * q;
    ret = q * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) /
      (((b[3] * r + b[2]) * r + b[1]) * r + one);
    return ret;;
  }
  else
  {
    if (q < zero)
      r = p;
    else
      r = one - p;

    if (r <= zero)
      return zero;
    r = sqrt (-log (r));
    if (r <= split2)
    {
      r = r - const2;
      ret = (((c[3] * r + c[2]) * r + c[1]) * r + c[0]) /
	((d[2] * r + d[1]) * r + one);
    }
    else
    {
      r = r - split2;
      ret = (((e[3] * r + e[2]) * r + e[1]) * r + e[0]) /
	((f[2] * r + f[1]) * r + one);
    }

    if (q < zero)
      ret = -ret;
    return ret;;
  }
}

double ppnd16 (double p)

/*-
 * algorithm as241  appl. statist. (1988) 37(3):
 *
 * produces the normal deviate z corresponding to a given lower
 * tail area of p; z is accurate to about 1 part in 10**16.
 *
 * the hash sums below are the sums of the mantissas of the
 * coefficients.   they are included for use in checking
 * transcription.
 */
{
  static double zero = 0.0, one = 1.0, half = 0.5;
  static double split1 = 0.425, split2 = 5.0;
  static double const1 = 0.180625, const2 = 1.6;

  /* coefficients for p close to 0.5 */
  static double a[8] = {
    3.3871328727963666080e0,
    1.3314166789178437745e+2,
    1.9715909503065514427e+3,
    1.3731693765509461125e+4,
    4.5921953931549871457e+4,
    6.7265770927008700853e+4,
    3.3430575583588128105e+4,
  2.5090809287301226727e+3};
  static double b[8] = {0.0,
    4.2313330701600911252e+1,
    6.8718700749205790830e+2,
    5.3941960214247511077e+3,
    2.1213794301586595867e+4,
    3.9307895800092710610e+4,
    2.8729085735721942674e+4,
  5.2264952788528545610e+3};

  /* hash sum ab    55.8831928806149014439 */
  /* coefficients for p not close to 0, 0.5 or 1. */
  static double c[8] = {
    1.42343711074968357734e0,
    4.63033784615654529590e0,
    5.76949722146069140550e0,
    3.64784832476320460504e0,
    1.27045825245236838258e0,
    2.41780725177450611770e-1,
    2.27238449892691845833e-2,
  7.74545014278341407640e-4};
  static double d[8] = {0.0,
    2.05319162663775882187e0,
    1.67638483018380384940e0,
    6.89767334985100004550e-1,
    1.48103976427480074590e-1,
    1.51986665636164571966e-2,
    5.47593808499534494600e-4,
  1.05075007164441684324e-9};

  /* hash sum cd    49.33206503301610289036 */
  /* coefficients for p near 0 or 1. */
  static double e[8] = {
    6.65790464350110377720e0,
    5.46378491116411436990e0,
    1.78482653991729133580e0,
    2.96560571828504891230e-1,
    2.65321895265761230930e-2,
    1.24266094738807843860e-3,
    2.71155556874348757815e-5,
  2.01033439929228813265e-7};
  static double f[8] = {0.0,
    5.99832206555887937690e-1,
    1.36929880922735805310e-1,
    1.48753612908506148525e-2,
    7.86869131145613259100e-4,
    1.84631831751005468180e-5,
    1.42151175831644588870e-7,
  2.04426310338993978564e-15};

  /* hash sum ef    47.52583317549289671629 */
  double q, r, ret;

  q = p - half;
  if (fabs (q) <= split1)
  {
    r = const1 - q * q;
    ret = q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3])
		 * r + a[2]) * r + a[1]) * r + a[0]) /
      (((((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3])
	 * r + b[2]) * r + b[1]) * r + one);
    return ret;
  }
  else
  {
    if (q < zero)
      r = p;
    else
      r = one - p;

    if (r <= zero)
      return zero;

    r = sqrt (-log (r));
    if (r <= split2)
    {
      r -= const2;
      ret = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3])
	       * r + c[2]) * r + c[1]) * r + c[0]) /
	(((((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3])
	   * r + d[2]) * r + d[1]) * r + one);
    }
    else
    {
      r -= split2;
      ret = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3])
	       * r + e[2]) * r + e[1]) * r + e[0]) /
	(((((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3])
	   * r + f[2]) * r + f[1]) * r + one);
    }

    if (q < zero)
      ret = -ret;
    return ret;
  }
}
/*-Algorithm AS 66
 * The Normal Integral, by I. D. Hill, 1973.
 * Applied Statistics 22(3):424-427.
 *
 * Translation to C by James Darrell McCauley, mccauley@ecn.purdue.edu.
 *
 * Calculates the upper or lower tail area of the standardized normal
 * curve corresponding to any given argument.
 *
 * x - the argument value
 * upper:  1 -> the area from x to \infty
 *         0 -> the area from -\infty to x
 *
 * Notes:
 * The constant LTONE should be set to the value at which the
 * lower tail area becomes 1.0 to the accuracy of the machine.
 * LTONE=(n+9)/3 gives the required value accurately enough, for a
 * machine that produces n decimal digits in its real numbers.
 *
 * The constant UTZERO should be set to the value at which the upper
 * tail area becomes 0.0 to the accuracy of the machine. This may be
 * taken as the value such that exp(-0.5 * UTZERO * UTZERO) /
 * (UTZERO * sqrt(2*PI)) is just greater than the smallest allowable
 * real numbers.
 */

#include<math.h>
#define LTONE 7.0
#define UTZERO 18.66

double alnorm (double x,int  upper)
{
  double ret, z, y;
  int up;

  up = upper;
  z = x;

  if (x < 0.0)
  {
    up = up == 0 ? 1 : 0;
    z = -x;
  }

  if (!(z <= LTONE || (up == 1 && z <= UTZERO)))
    ret = 0.0;
  else
  {
    y = 0.5 * z * z;
    if (z <= 1.28)
      ret = 0.5 - z * (0.398942280444 - 0.399903438504 * y /
		       (y + 5.75885480458 - 29.8213557808 /
			(y + 2.62433121679 + 48.6959930692 /
			 (y + 5.92885724438))));
    else
      ret = 0.398942280385 * exp (-y) /
	(z - 3.8052e-8 + 1.00000615302 /
	 (z + 3.98064794e-4 + 1.98615381364 /
	  (z - 0.151679116635 + 5.29330324926 /
	   (z + 4.8385912808 - 15.1508972451 /
	    (z + 0.742380924027 + 30.789933034 /
	     (z + 3.99019417011))))));
  }
  if (up == 0)
    ret = 1.0 - ret;
  return ret;
}

double *chi_square_exp (double* x,int n)
{
  static double y[2];
  double mean = 0.0, sum3 = 0.0, *v;
  int i, j, k, *f;

  k = (int)(4.0 * pow (0.75 * (n - 1.0) * (n - 1.0), 0.2));//bbcrevisit

  while ((double) (n/k) < 5.0)
   --k;

  if ((f = (int *) calloc (k, sizeof (int))) == NULL)
    fprintf (stderr, "Memory error in chi_square\n"), exit (-1);
  if ((v = (double *) malloc ( (k+1) * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in chi_square\n"), exit (-1);

  for (i = 0; i < n; ++i)
    mean += x[i];
  mean = n/mean;

  v[0]=0.0;
  for (i = 1; i < k; ++i)
    v[i] = - log (1.0 - (double) i / k) / mean;
  v[k]=1e9;

  for (i = 0; i < n; ++i)
  {
    j=0;
    while(j<k)
    {
      if (x[i] > v[j] && x[i] <= v[j+1])
      {
	f[j]++;
        j=k;
      }
     j++;
    }
  }

  for (i = 0; i < k; ++i)
    sum3 += f[i] * f[i];

  y[0] = sum3 * k / n - n;
  y[1] = (double) k - 2.0;
#ifdef NOISY
  printf ("  TEST21 CS(E)  =%10.4f   DOF    =%10.4f\n", y[0], y[1]);
#endif /* NOISY */
  free(f);
  free(v);
  return y;
}

double *chi_square (double* x,int  n)
{
  static double y[2];
  double mean = 0.0, sdx = 0.0, sum3 = 0.0, *v;
  int i, j, k, *f;

  k = (int)(4.0 * pow (0.75 * (n - 1.0) * (n - 1.0), 0.2));

  while ((double) (n/k) < 5.0)
   --k;

  if ((f = (int *) calloc (k, sizeof (int))) == NULL)
    fprintf (stderr, "Memory error in chi_square\n"), exit (-1);
  if ((v = (double *) malloc ( (k+1) * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in chi_square\n"), exit (-1);

  for (i = 0; i < n; ++i)
  {
    mean += x[i];
    sdx += x[i] * x[i];
  }
  sdx = sqrt ((n * sdx - mean * mean) / (n * (n - 1.0)));
  mean /= n;

  v[0]=-1e9;
  for (i = 1; i < k; ++i)
    v[i] = mean + xinormal ((double) i / k) * sdx;
  v[k]=1e9;

  for (i = 0; i < n; ++i)
  {
    j=0;
    while(j<k)
    {
      if (x[i] > v[j] && x[i] <= v[j+1])
      {
	f[j]++;
        j=k;
      }
     j++;
    }
  }

  for (i = 0; i < k; ++i)
    sum3 += f[i] * f[i];

  y[0] = sum3 * k / n - n;
  y[1] = (double) k - 3.0;
#ifdef NOISY
  printf ("  TEST12 CS(N)  =%10.4f   DOF    =%10.4f\n", y[0], y[1]);
#endif /* NOISY */
  free(f);
  free(v);
  return y;
}

double *cramer_von_mises_exp (double* x,int  n)
{
  static double y[2];
  double *xcopy, mean=0.0, fx, fn2, sum4 = 0.0;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in cramer_von_mises_exp\n"), exit (-1);

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
  }
  mean /= n;

  qsort (xcopy, n, sizeof (double), dcmp);

  for (i = 0; i < n; ++i)
  {
    /*-
    a = (2 * i + 1) * log (fx);
    b = (2 * i + 1) * (xcopy[n-i-1] * (-1.0 / mean));
    sum3 += a + b;
    */
    fx = 1 - exp (xcopy[i] * (-1.0 / mean));
    fn2 = (double) (2.0 * i + 1) / (2 * n);
    sum4 += (fx - fn2) * (fx - fn2);
  }

  /*-
  cvm = 1.0 / (n * 12) + sum4;
  cvmod = cvm * (0.16 / n + 1.0);
  */
  y[0] = (1.0 / (n * 12) + sum4)*(0.16 / n + 1.0);

#ifdef NOISY
  printf ("  TEST16 CVM(E) =%10.4f\n", y[0]);
#endif				/* NOISY */

  free (xcopy);
  return y;
}

double *cramer_von_mises (double * x,int  n)
{
  int i ;
  static double y[2];
  double mean = 0.0, sdx = 0.0, fx, sqrt2, *xcopy;

  sqrt2 = sqrt ((double) 2.0);
  y[1] = 0.0;
  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in cramer_von_mises\n"), exit (-1);

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
    sdx += x[i] * x[i];
  }
  sdx = sqrt ((n * sdx - mean * mean) / (n * (n - 1.0)));
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    fx = 0.5 + normp ((xcopy[i] - mean) / sdx / sqrt2) / 2.0;
    if (fx <= 1e-5)
      fx = 1e-5;
    if (fx >= 0.99999)
      fx = 0.99999;
    fx -= (2.0 * i + 1.0) / (2.0 * n);
    y[1] += fx * fx;
  }
  y[1] += 1.0 / (double) (n * 12);
  y[0] = y[1] * (0.5 / n + 1.0);
#ifdef NOISY
  printf ("  TEST9  CVM(N) =%10.4f\n", y[0]);
#endif				/* NOISY */
  free(xcopy);
  return y;
}
//bbcrevisit - not enough cmonguys
#define PI 3.141592654

double *dagostino_d (double* x,int  n)
{
  int i;
  static double y[2];
  double d, s, t = 0., *xcopy, m2, s1 = 0., s2, mn = 0.0;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
  {
    fprintf (stderr, "Memory allocation error\n");
    exit (-1);
  }
  for (i = 0; i < n; ++i)
    xcopy[i] = x[i];
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    t += xcopy[i] * ((i + 1) - 0.5 * (n + 1));
    mn += xcopy[i];
  }
  m2 = mn / n;
  for (i = 0; i < n; ++i)
    s1 += (xcopy[i] - m2) * (xcopy[i] - m2);
  s2 = s1 / n;
  s = sqrt (s2);
  d = t / (n * n * s);
  /* y[0] = (d - 1. / (2*sqrt (PI))) * sqrt ((double)n) / 0.02998598; */
  y[0] = d;
  y[1] = sqrt((double)n)*(y[0]-0.28209479)/0.02998598;

#ifdef NOISY
  printf ("  TEST4  DAGN   =%10.4f\n", y[0]);
#endif /* NOISY */
  return y;
}
/* this if the comparison function for the qsort */
int dblcomp (double* i,double* j)
  
{
  if (*i - *j < 0)
    return -1;
  else if (*i - *j > 0)
    return 1;
  else
    return 0;
}

double *dmax (double *x,int n)

{
  static double y[2];
  double *xcopy, sqrt2, sqrtn, mean = 0.0, sdx = 0.0, fx;
  double dp, dp_max, dm, dm_max;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in dmax\n"), exit (-1);

  sqrt2 = sqrt ((double) 2.0);
  sqrtn = sqrt ((double) n);

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
    sdx += x[i] * x[i];
  }
  sdx = sqrt ((n * sdx - mean * mean) / (n * (n - 1.0)));
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = (xcopy[i] - mean) / sdx;
    fx = 0.5 + normp (xcopy[i] / sqrt2) / 2.0;
    if (fx <= 1e-5)
      fx = 1e-5;
    if (fx >= 0.99999)
      fx = 0.99999;
    dp = (double) (i + 1) / (double) n - fx;
    dm = fx - i / (double) n;
    if (i == 0 || dp > dp_max)
      dp_max = dp;
    if (i == 0 || dm > dm_max)
      dm_max = dm;
  }
  y[0] = dp_max;
  y[1] = dm_max;
  free (xcopy);
  return y;
}

double *dmax_exp (double *x, int n)
{
  static double y[2];
  double mean=0.0, zmax, tmax, *xcopy, t, z, fx;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in dmax_exp\n"), exit (-1);
 
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
  }
  mean /=n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    fx = 1-exp(-xcopy[i]/mean);
    z = (double) (i+1) / (double) n - fx;
    t = fx - (double) i / (double) n;
    if (i == 0 || z > zmax)
      zmax = z;
    if (i == 0 || t > tmax)
      tmax = t;
  }
  y[0]=zmax;
  y[1]=tmax;
  free(xcopy);
  return y;
}

double *durbins_exact ( double *x,  int n)
{
  static double y[2];
  double *xcopy, sumx = 0.0, sumx2 = 0.0, s2, *b, *c, *g, *z, sqrt2;
 
  int i, j;

  if ((b = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in durbins_exact\n"), exit (-1);
  if ((c = (double *) malloc ((n + 1) * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in durbins_exact\n"), exit (-1);
  if ((g = (double *) malloc ((n + 1) * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in durbins_exact\n"), exit (-1);
  if ((z = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in durbins_exact\n"), exit (-1);
  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in durbins_exact\n"), exit (-1);

  sqrt2 = sqrt ((double) 2.0);
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    sumx += x[i];
    sumx2 += x[i] * x[i];
  }
  s2 = sqrt ((sumx2 - sumx * sumx / n) / (n - 1));
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = (xcopy[i] - sumx / n) / s2;
    b[i] = 0.5 + normp (xcopy[i] / sqrt2) / 2.0;
  }
  qsort (b, n, sizeof (double), dcmp);
  for (i = 1; i < n; ++i)
    c[i] = b[i] - b[i - 1];
  c[0] = b[0];
  c[n] = 1 - b[n - 1];
  qsort (c, n + 1, sizeof (double), dcmp);
  for (j = 1; j <= n; ++j)
    g[j] = (n + 1 - j) * (c[j] - c[j - 1]);
  g[0] = (n + 1) * c[0];
  g[n] = c[n] - c[n - 1];

  for (i = 0; i < n; ++i)
  {
    z[i] = 0.0;
    for (j = 0; j <= i; ++j)
      z[i] += g[j];
    z[i] = (i + 1.0) / n - z[i];
  }

  qsort (z, n, sizeof (double), dcmp);
  y[0] = z[n - 1];
  y[1] = sqrt((double)n)*z[n - 1];

#ifdef NOISY
  printf ("  TEST7  DRB(N) =%10.4f\n", y[0]);
#endif				/* NOISY */
  free (b);
  free (c);
  free (g);
  free (xcopy);
  free (z);
  return y;
}

double enormp (double x)
  
{
  double x1, x2, x3, x4, ret_val;
  static double xp[5] = {7.7105849500132e-5, -0.00133733772997339,
  0.0323076579225834, 0.0479137145607681, 0.128379167095513};
  static double xq[3] = {0.00301048631703895, 0.0538971687740286,
  0.375795757275549};
  static double xr[8] = {-1.36864857382717e-7, 0.564195517478974,
    7.21175825088309, 43.1622272220567, 152.98928504694,
  339.320816734344, 451.918953711873, 300.459261020162};
  static double xs[8] = {1.0, 12.7827273196294, 77.0001529352295,
    277.585444743988, 638.980264465631, 931.35409485061,
  790.950925327898, 300.459260956983};
  static double xt[5] = {2.10144126479064, 26.2370141675169,
  21.3688200555087, 4.6580782871847, 0.282094791773523};
  static double xu[4] = {94.153775055546, 187.11481179959,
  99.0191814623914, 18.0124575948747};
  double yy1, yy2;

  x3 = (double) 0.564189583547756;
  x1 = fabs (x);

  if (x1 <= .5)
  {
    x4 = x * x;
    yy1 = (((xp[0] * x4 + xp[1]) * x4 + xp[2]) * x4 + xp[3]) * x4 + xp[4] + 1.;
    yy2 = ((xq[0] * x4 + xq[1]) * x4 + xq[2]) * x4 + 1.;
    ret_val = x * (yy1 / yy2);
  }
  else if (x1 <= 4.)
  {
    yy1 = ((((((xr[0] * x1 + xr[1]) * x1 + xr[2]) * x1 + xr[3])
	     * x1 + xr[4]) * x1 + xr[5]) * x1 + xr[6]) * x1 + xr[7];
    yy2 = ((((((xs[0] * x1 + xs[1]) * x1 + xs[2]) * x1 + xs[3])
	     * x1 + xs[4]) * x1 + xs[5]) * x1 + xs[6]) * x1 + xs[7];
    ret_val = 1. - exp (-x * x) * yy1 / yy2;
    if (x < 0.)
      ret_val *= -1.0;
  }
  else
  {
    x2 = x * x;
    x4 *= 1.; /* huh, what? See original FORTRAN */
    x4 = 0.5;
    yy1 = (((xt[0] * x4 + xt[1]) * x4 + xt[2]) * x4 + xt[3]) * x4 + xt[4];
    yy2 = (((xu[0] * x4 + xu[1]) * x4 + xu[2]) * x4 + xu[3]) * x4 + 1.;
    ret_val = x3 / x1 - yy1 * x1 / (x2 * yy2);
    ret_val = 1. - exp (-x2) * ret_val;
    if (x < 0.)
      ret_val = -ret_val;
  }
  return ret_val;
}

double *extreme (  double * x,   int n)

{
  int i;
  static double y[2];
  double min, max, sum1 = 0.;

  min = max = x[0];
  for (i = 0; i < n; ++i)
  {
    sum1 += x[i];
    if (min > x[i])
      min = x[i];
    if (max < x[i])
      max = x[i];
  }
  sum1 /= n;
  y[0] = max - sum1;
  y[1] = min - sum1;

#ifdef NOISY
  printf ("  TEST3  U(N)   =%10.4f   U(1)   =%10.4f\n", y[0], y[1]);
#endif				/* NOISY */
  return y;
}

double *geary_test (double *x,  int  n)

{
  int i;
  static double y[2];
  double diff, s=0.0, mean = 0.0;

  y[0] = 0.0;
  for (i = 0; i < n; ++i)
    mean += x[i];
  mean /= n;

  for (i = 0; i < n; ++i)
  {
    diff = x[i] - mean;
    y[0] += fabs (diff);
    s += diff * diff;
  }

  s *= n;
  y[0] /= sqrt (s);
  y[1] = (y[0] - 0.7979) * sqrt ((double) n) / 0.2123;

#ifdef NOISY
  printf ("  TEST2  GTN    =%10.4f   Z(GTN) =%10.4f\n", y[0], y[1]);
#endif				/* NOISY */
  return y;
}

double *kotz_families (double *x,int n)

{
  static double y[2];
  int i;
  double a1, b1, a2, b3, c1, c2, c3, c4, c5, c6, lx;
  double sum1 = 0.0, sum2 = 0.0, sum4 = 0.0;

  for (i = 0; i < n; ++i)
  {
    sum1 += x[i];
    sum2 += log (x[i]);
  }
  b1 = sum1 / n;
  a1 = sum2 / n;
  for (i = 0; i < n; ++i)
  {
    lx = log (x[i]);
    sum4 += (lx - a1) * (lx - a1);
  }
  a2 = sum4 / n;
  b3 = exp (a1 * 2 + a2) * (exp (a2) - 1);
  c1 = log (a2 / b3);
  c2 = (exp (a2 * 4) + exp (a2 * 3) * 2 - 4) / 4 - a2 + exp (a2) * 0.75;
  c3 = a2 * (exp (a2) * 2 - 1) * (exp (a2) * 2 - 1);
  c4 = (exp (a2) - 1) * 2 * (exp (a2) - 1);
  c5 = c3 / c4;
  if (c2 < c5)
  {
#ifdef NOISY
    printf ("  WARNING!!! STATISTICS FOR THE NEXT TEST WILL\n");
    printf ("  NOT BE CALCULATED DUE TO SMALL LOGVARIANCE\n");
#endif				/* NOISY */
    y[0] = 999999999.;
  }
  else
  {
    c6 = sqrt (c2 - c5) * 2.0 * sqrt ((double) n);
    y[0] = c1 / c6;
  }
#ifdef NOISY
  printf ("  TEST24 KT(LN) =%10.4f\n", y[0]);
#endif				/* NOISY */
  return y;
}

double *kolmogorov_smirnov_exp (double* x,int  n)
{
  static double y[2];
  double *d, sqrtn;

  d=dmax_exp(x,n);
  sqrtn = sqrt ((double) n);
  y[1] = (d[0] > d[1]) ? d[0] : d[1];
  y[0] = (y[1] - 0.2 / n) * (sqrtn + 0.5 / sqrtn + 0.26);
#ifdef NOISY
  printf ("  TEST17 KSD(E) =%10.4f\n", y[0]);
#endif				/* NOISY */
  return y;
}

double *kolmogorov_smirnov (double *x,int n)

{
  static double y[2];
  double *d,sqrtn;

  sqrtn = sqrt ((double) n);
  d=dmax(x,n);

  y[1] = (d[0] > d[1]) ? d[0] : d[1];
  y[0] = y[1] * (sqrtn + 0.85 / sqrtn - 0.01);
#ifdef NOISY
  printf ("  TEST10 KSD(N) =%10.4f\n", y[0]);
  printf ("  TEST11 KSD    =%10.4f\n", y[1]);
#endif				/* NOISY */
  return y;
}

double *kuipers_v_exp (  double *x ,  int n)
{
  static double y[2];
  double *d,r;

  d = dmax_exp (x, n);
  r = sqrt ((double) n);
  y[1] = d[0] + d[1];
  y[0] = (y[1] - 0.2 / n) * (r + 0.35 / r + 0.24);
#ifdef NOISY
  printf ("  TEST18 KV(E)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  return y;
}
double *kuipers_v (double * x,int n)

{
  static double y[2]; 
  double *d, sqrtn;
  
  sqrtn=sqrt((double)n);
  d=dmax(x,n);

  y[1] = d[0] + d[1];
  y[0] = y[1] * (sqrtn + 0.05 + 0.82 / sqrtn);
#ifdef NOISY
  printf ("  TEST5  KV(N)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  return y;
}

/*-
 *     SUBROUTINE NORMP(Z, P, Q, PDF)
 *
 *     Normal distribution probabilities accurate to 1.e-15.
 *     Z = no. of standard deviations from the mean.
 *     P, Q = probabilities to the left & right of Z.   P + Q = 1.
 *     PDF = the probability density.
 *
 *     Based upon algorithm 5666 for the error function, from:
 *     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
 *
 *     Programmer: Alan Miller
 *
 *     Latest revision - 30 March 1986
 *
 */

/* Conversion to C by James Darrell McCauley, 24 September 1994 */

double normp (double z)
{
  static double p[7] = {220.2068679123761, 221.2135961699311,
    112.079291497870, 33.91286607838300, 6.37396220353165,
  0.7003830644436881, 0.352624965998910e-1};
  static double q[8] = {440.4137358247522, 793.8265125199484,
    637.3336333788311, 296.5642487796737, 86.78073220294608,
  16.06417757920695, 1.755667163182642, 0.8838834764831844e-1};
  static double cutoff = 7.071, root2pi = 2.506628274631001;
  double zabs, expntl;
  double pee, queue, pdf;

  zabs = fabs (z);

  if (zabs > 37.0)
  {
    pdf = 0.0;
    if (z > 0.0)
    {
      pee = 1.0;
      queue = 0.0;
    }
    else
    {
      pee = 0.0;
      queue = 1.0;
    }
    return pee;
  }

  expntl = exp (-0.5 * zabs * zabs);
  pdf = expntl / root2pi;

    if (zabs < cutoff)
    pee = expntl * ((((((p[6] * zabs + p[5]) * zabs + p[4])
	      * zabs + p[3]) * zabs + p[2]) * zabs + p[1]) * zabs + p[0])
      / (((((((q[7] * zabs + q[6]) * zabs + q[5]) * zabs + q[4])
	    * zabs + q[3]) * zabs + q[2]) * zabs + q[1]) * zabs + q[0]);
  else
    pee = pdf / (zabs + 1.0 / (zabs + 2.0 / (zabs + 3.0 / (zabs + 4.0 /
						       (zabs + 0.65)))));

  if (z < 0.0)
    queue = 1.0 - pee;
  else
  {
    queue = pee;
      pee = 1.0 - queue;
  }
  return pee;
}

double *omnibus_moments (double * x, int n)
{
  double diff, mean = 0., fssm, tssm, sum_cube = 0., sum_four = 0.,
   sum_sq = 0.;
  static double y[2];
  int i;

  for (i = 0; i < n; ++i)
    mean += x[i];
  mean /= n;

  for (i = 0; i < n; ++i)
  {
    diff = x[i] - mean;
    sum_sq += diff * diff;
    sum_cube += diff * diff * diff;
    sum_four += diff * diff * diff * diff;
  }

/*
printf("n %d x-bar %g sum^2 %g sum^3 %g sum^4 %g \n",n,mean,sum_sq,sum_cube,sum_four);
*/
  tssm = sqrt ((double) n) * sum_cube / pow (sum_sq, 1.5);
  fssm = n * sum_four / (sum_sq * sum_sq);

#ifdef NOISY
  printf ("          TESTS OF COMPOSITE DISTRIBUTIONAL HYPOTHESES\n");
  printf ("  TEST1  TSM    =%10.4f   FSM    =%10.4f\n", tssm, fssm);
#endif				/* NOISY */

  y[0] = tssm;
  y[1] = fssm;
  return y;
}

/*-
 * driver program for AS 181: Royston's extension of the Shapiro-Wilk
 * W statistic to n=2000
 * needs as181.c as177.c as241.c dcmp.c as66.c 
 */

double *royston(double * x,int n)

{  static double y[2];
  double *a, eps, w, pw, mean=0, ssq=0, *xcopy;
  int i,  n2; 
int* ifault=new int;
*ifault=0;

  n2=(int)floor((double)n/2);

#ifndef lint
  if ((a=(double *) malloc(n2*sizeof(double)))==NULL)
    fprintf (stderr, "Memory error in royston\n"), exit (-1);
  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in shapiro_wilk\n"), exit (-1);
#endif // lint 
 

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
  }
  mean /=n ;

  qsort (xcopy, n, sizeof (double), dcmp);

  for (i = 0; i < n; ++i)
    ssq += (mean-x[i])*(mean-x[i]);

  wcoef (a, n, n2, &eps, ifault);

  if (*ifault==0)
    wext (xcopy, n, ssq, a, n2, eps, &w, &pw, ifault);
  else
  {
    fprintf (stderr, "Error in wcoef()\n");
    return (double *) NULL;
  }

  if (ifault==0)
  {
    y[0]=w;
    y[1]=pw;
  }
  else
  {
    fprintf (stderr, "Error in wcoef()\n");
    return (double *) NULL;
  }

  free(a);
  free(xcopy);
  return y;
  
}

double *shapiro_wilk (double * x,  int  n)
{
  static double y[2];
  double a[25], s2, *xcopy;
  double sumb = 0.0, sumx = 0.0, sumx2 = 0.0;
  int i, k;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in shapiro_wilk\n"), exit (-1);

  k = n / 2;
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    sumx += x[i];
    sumx2 += x[i] * x[i];
  }
  s2 = sumx2 - sumx * sumx / n;
  qsort (xcopy, n, sizeof (double), dcmp);

  if (n == 3)
  {
    a[0] = (double) .7071;
  }
  else if (n == 4)
  {
    a[0] = (double) .6872;
    a[1] = (double) .1677;
  }
  else if (n == 5)
  {
    a[0] = (double) .6646;
    a[1] = (double) .2413;
  }
  else if (n == 6)
  {
    a[0] = (double) .6431;
    a[1] = (double) .2806;
    a[2] = (double) .0875;
  }
  else if (n == 7)
  {
    a[0] = (double) .6233;
    a[1] = (double) .3031;
    a[2] = (double) .1401;
  }
  else if (n == 8)
  {
    a[0] = (double) .6052;
    a[1] = (double) .3164;
    a[2] = (double) .1743;
    a[3] = (double) .0561;
  }
  else if (n == 9)
  {
    a[0] = (double) .5888;
    a[1] = (double) .3244;
    a[2] = (double) .1976;
    a[3] = (double) .0947;
  }
  else if (n == 10)
  {
    a[0] = (double) .5739;
    a[1] = (double) .3291;
    a[2] = (double) .2141;
    a[3] = (double) .1224;
    a[4] = (double) .0399;
  }
  else if (n == 11)
  {
    a[0] = (double) .5601;
    a[1] = (double) .3315;
    a[2] = (double) .226;
    a[3] = (double) .1429;
    a[4] = (double) .0695;
  }
  else if (n == 12)
  {
    a[0] = (double) .5475;
    a[1] = (double) .3325;
    a[2] = (double) .2347;
    a[3] = (double) .1586;
    a[4] = (double) .0922;
    a[5] = (double) .0303;
  }
  else if (n == 13)
  {
    a[0] = (double) .5359;
    a[1] = (double) .3325;
    a[2] = (double) .2412;
    a[3] = (double) .1707;
    a[4] = (double) .1099;
    a[5] = (double) .0539;
  }
  else if (n == 14)
  {
    a[0] = (double) .5251;
    a[1] = (double) .3318;
    a[2] = (double) .246;
    a[3] = (double) .1802;
    a[4] = (double) .124;
    a[5] = (double) .0727;
    a[6] = (double) .024;
  }
  else if (n == 15)
  {
    a[0] = (double) .515;
    a[1] = (double) .3306;
    a[2] = (double) .2495;
    a[3] = (double) .1878;
    a[4] = (double) .1353;
    a[5] = (double) .088;
    a[6] = (double) .0433;
  }
  else if (n == 16)
  {
    a[0] = (double) .5056;
    a[1] = (double) .329;
    a[2] = (double) .2521;
    a[3] = (double) .1939;
    a[4] = (double) .1447;
    a[5] = (double) .1005;
    a[6] = (double) .0593;
    a[7] = (double) .0196;
  }
  else if (n == 17)
  {
    a[0] = (double) .4968;
    a[1] = (double) .3273;
    a[2] = (double) .254;
    a[3] = (double) .1988;
    a[4] = (double) .1524;
    a[5] = (double) .1109;
    a[6] = (double) .0725;
    a[7] = (double) .0359;
  }
  else if (n == 18)
  {
    a[0] = (double) .4886;
    a[1] = (double) .3253;
    a[2] = (double) .2553;
    a[3] = (double) .2027;
    a[4] = (double) .1587;
    a[5] = (double) .1197;
    a[6] = (double) .0837;
    a[7] = (double) .0496;
    a[8] = (double) .0163;
  }
  else if (n == 19)
  {
    a[0] = (double) .4808;
    a[1] = (double) .3232;
    a[2] = (double) .2561;
    a[3] = (double) .2059;
    a[4] = (double) .1641;
    a[5] = (double) .1271;
    a[6] = (double) .0932;
    a[7] = (double) .0612;
    a[8] = (double) .0303;
  }
  else if (n == 20)
  {
    a[0] = (double) .4734;
    a[1] = (double) .3211;
    a[2] = (double) .2565;
    a[3] = (double) .2085;
    a[4] = (double) .1686;
    a[5] = (double) .1334;
    a[6] = (double) .1013;
    a[7] = (double) .0711;
    a[8] = (double) .0422;
    a[9] = (double) .014;
  }
  else if (n == 21)
  {
    a[0] = (double) .4643;
    a[1] = (double) .3185;
    a[2] = (double) .2578;
    a[3] = (double) .2119;
    a[4] = (double) .1736;
    a[5] = (double) .1399;
    a[6] = (double) .1092;
    a[7] = (double) .0804;
    a[8] = (double) .053;
    a[9] = (double) .0263;
  }
  else if (n == 22)
  {
    a[0] = (double) .459;
    a[1] = (double) .3156;
    a[2] = (double) .2571;
    a[3] = (double) .2131;
    a[4] = (double) .1764;
    a[5] = (double) .1443;
    a[6] = (double) .115;
    a[7] = (double) .0878;
    a[8] = (double) .0618;
    a[9] = (double) .0368;
    a[10] = (double) .0122;
  }
  else if (n == 23)
  {
    a[0] = (double) .4542;
    a[1] = (double) .3126;
    a[2] = (double) .2563;
    a[3] = (double) .2139;
    a[4] = (double) .1787;
    a[5] = (double) .148;
    a[6] = (double) .1201;
    a[7] = (double) .0941;
    a[8] = (double) .0696;
    a[9] = (double) .0459;
    a[10] = (double) .0228;
  }
  else if (n == 24)
  {
    a[0] = (double) .4493;
    a[1] = (double) .3098;
    a[2] = (double) .2554;
    a[3] = (double) .2145;
    a[4] = (double) .1807;
    a[5] = (double) .1512;
    a[6] = (double) .1245;
    a[7] = (double) .0997;
    a[8] = (double) .0764;
    a[9] = (double) .0539;
    a[10] = (double) .0321;
    a[11] = (double) .0107;
  }
  else if (n == 25)
  {
    a[0] = (double) .445;
    a[1] = (double) .3069;
    a[2] = (double) .2543;
    a[3] = (double) .2148;
    a[4] = (double) .1822;
    a[5] = (double) .1539;
    a[6] = (double) .1283;
    a[7] = (double) .1046;
    a[8] = (double) .0823;
    a[9] = (double) .061;
    a[10] = (double) .0403;
    a[11] = (double) .02;
  }
  else if (n == 26)
  {
    a[0] = (double) .4407;
    a[1] = (double) .3043;
    a[2] = (double) .2533;
    a[3] = (double) .2151;
    a[4] = (double) .1836;
    a[5] = (double) .1563;
    a[6] = (double) .1316;
    a[7] = (double) .1089;
    a[8] = (double) .0876;
    a[9] = (double) .0672;
    a[10] = (double) .0476;
    a[11] = (double) .0284;
    a[12] = (double) .0094;
  }
  else if (n == 27)
  {
    a[0] = (double) .4366;
    a[1] = (double) .3018;
    a[2] = (double) .2522;
    a[3] = (double) .2152;
    a[4] = (double) .1848;
    a[5] = (double) .1584;
    a[6] = (double) .1346;
    a[7] = (double) .1128;
    a[8] = (double) .0923;
    a[9] = (double) .0728;
    a[10] = (double) .054;
    a[11] = (double) .0358;
    a[12] = (double) .0178;
  }
  else if (n == 28)
  {
    a[0] = (double) .4328;
    a[1] = (double) .2992;
    a[2] = (double) .251;
    a[3] = (double) .2151;
    a[4] = (double) .1857;
    a[5] = (double) .1601;
    a[6] = (double) .1372;
    a[7] = (double) .1162;
    a[8] = (double) .0965;
    a[9] = (double) .0778;
    a[10] = (double) .0598;
    a[11] = (double) .0424;
    a[12] = (double) .0253;
    a[13] = (double) .0084;
  }
  else if (n == 29)
  {
    a[0] = (double) .4291;
    a[1] = (double) .2968;
    a[2] = (double) .2499;
    a[3] = (double) .215;
    a[4] = (double) .1864;
    a[5] = (double) .1616;
    a[6] = (double) .1395;
    a[7] = (double) .1192;
    a[8] = (double) .1002;
    a[9] = (double) .0822;
    a[10] = (double) .065;
    a[11] = (double) .0483;
    a[12] = (double) .032;
    a[13] = (double) .0159;
  }
  else if (n == 30)
  {
    a[0] = (double) .4254;
    a[1] = (double) .2944;
    a[2] = (double) .2487;
    a[3] = (double) .2148;
    a[4] = (double) .187;
    a[5] = (double) .163;
    a[6] = (double) .1415;
    a[7] = (double) .1219;
    a[8] = (double) .1036;
    a[9] = (double) .0862;
    a[10] = (double) .0697;
    a[11] = (double) .0537;
    a[12] = (double) .0381;
    a[13] = (double) .0227;
    a[14] = (double) .0076;
  }
  else if (n == 31)
  {
    a[0] = (double) .422;
    a[1] = (double) .2921;
    a[2] = (double) .2475;
    a[3] = (double) .2145;
    a[4] = (double) .1874;
    a[5] = (double) .1641;
    a[6] = (double) .1433;
    a[7] = (double) .1243;
    a[8] = (double) .1066;
    a[9] = (double) .0899;
    a[10] = (double) .0739;
    a[11] = (double) .0585;
    a[12] = (double) .0435;
    a[13] = (double) .0289;
    a[14] = (double) .0144;
  }
  else if (n == 32)
  {
    a[0] = (double) .4188;
    a[1] = (double) .2898;
    a[2] = (double) .2463;
    a[3] = (double) .2141;
    a[4] = (double) .1878;
    a[5] = (double) .1651;
    a[6] = (double) .1449;
    a[7] = (double) .1265;
    a[8] = (double) .1093;
    a[9] = (double) .0931;
    a[10] = (double) .0777;
    a[11] = (double) .0629;
    a[12] = (double) .0485;
    a[13] = (double) .0344;
    a[14] = (double) .0206;
    a[15] = (double) .0068;
  }
  else if (n == 33)
  {
    a[0] = (double) .4156;
    a[1] = (double) .2876;
    a[2] = (double) .2451;
    a[3] = (double) .2137;
    a[4] = (double) .188;
    a[5] = (double) .166;
    a[6] = (double) .1463;
    a[7] = (double) .1284;
    a[8] = (double) .1118;
    a[9] = (double) .0961;
    a[10] = (double) .0812;
    a[11] = (double) .0669;
    a[12] = (double) .053;
    a[13] = (double) .0395;
    a[14] = (double) .0262;
    a[15] = (double) .0131;
  }
  else if (n == 34)
  {
    a[0] = (double) .4127;
    a[1] = (double) .2854;
    a[2] = (double) .2439;
    a[3] = (double) .2132;
    a[4] = (double) .1882;
    a[5] = (double) .1667;
    a[6] = (double) .1475;
    a[7] = (double) .1301;
    a[8] = (double) .114;
    a[9] = (double) .0988;
    a[10] = (double) .0844;
    a[11] = (double) .0706;
    a[12] = (double) .0572;
    a[13] = (double) .0441;
    a[14] = (double) .0314;
    a[15] = (double) .0187;
    a[16] = (double) .0062;
  }
  else if (n == 35)
  {
    a[0] = (double) .4096;
    a[1] = (double) .2834;
    a[2] = (double) .2427;
    a[3] = (double) .2127;
    a[4] = (double) .1883;
    a[5] = (double) .1673;
    a[6] = (double) .1487;
    a[7] = (double) .1317;
    a[8] = (double) .116;
    a[9] = (double) .1013;
    a[10] = (double) .0873;
    a[11] = (double) .0739;
    a[12] = (double) .061;
    a[13] = (double) .0484;
    a[14] = (double) .0361;
    a[15] = (double) .0239;
    a[16] = (double) .0119;
  }
  else if (n == 36)
  {
    a[0] = (double) .4068;
    a[1] = (double) .2813;
    a[2] = (double) .2415;
    a[3] = (double) .2121;
    a[4] = (double) .1883;
    a[5] = (double) .1678;
    a[6] = (double) .1496;
    a[7] = (double) .1331;
    a[8] = (double) .1179;
    a[9] = (double) .1036;
    a[10] = (double) .09;
    a[11] = (double) .077;
    a[12] = (double) .0645;
    a[13] = (double) .0523;
    a[14] = (double) .0404;
    a[15] = (double) .0287;
    a[16] = (double) .0172;
    a[17] = (double) .0057;
  }
  else if (n == 37)
  {
    a[0] = (double) .404;
    a[1] = (double) .2794;
    a[2] = (double) .2403;
    a[3] = (double) .2116;
    a[4] = (double) .1883;
    a[5] = (double) .1683;
    a[6] = (double) .1505;
    a[7] = (double) .1344;
    a[8] = (double) .1196;
    a[9] = (double) .1056;
    a[10] = (double) .0924;
    a[11] = (double) .0798;
    a[12] = (double) .0677;
    a[13] = (double) .0559;
    a[14] = (double) .0444;
    a[15] = (double) .0331;
    a[16] = (double) .022;
    a[17] = (double) .011;
  }
  else if (n == 38)
  {
    a[0] = (double) .4015;
    a[1] = (double) .2774;
    a[2] = (double) .2391;
    a[3] = (double) .211;
    a[4] = (double) .1881;
    a[5] = (double) .1686;
    a[6] = (double) .1513;
    a[7] = (double) .1356;
    a[8] = (double) .1211;
    a[9] = (double) .1075;
    a[10] = (double) .0947;
    a[11] = (double) .0824;
    a[12] = (double) .0706;
    a[13] = (double) .0592;
    a[14] = (double) .0481;
    a[15] = (double) .0372;
    a[16] = (double) .0264;
    a[17] = (double) .0158;
    a[18] = (double) .0053;
  }
  else if (n == 39)
  {
    a[0] = (double) .3989;
    a[1] = (double) .2755;
    a[2] = (double) .238;
    a[3] = (double) .2104;
    a[4] = (double) .188;
    a[5] = (double) .1689;
    a[6] = (double) .152;
    a[7] = (double) .1366;
    a[8] = (double) .1225;
    a[9] = (double) .1092;
    a[10] = (double) .0967;
    a[11] = (double) .0848;
    a[12] = (double) .0733;
    a[13] = (double) .0622;
    a[14] = (double) .0515;
    a[15] = (double) .0409;
    a[16] = (double) .0305;
    a[17] = (double) .0203;
    a[18] = (double) .0101;
  }
  else if (n == 40)
  {
    a[0] = (double) .3964;
    a[1] = (double) .2737;
    a[2] = (double) .2368;
    a[3] = (double) .2098;
    a[4] = (double) .1878;
    a[5] = (double) .1691;
    a[6] = (double) .1526;
    a[7] = (double) .1376;
    a[8] = (double) .1237;
    a[9] = (double) .1108;
    a[10] = (double) .0986;
    a[11] = (double) .087;
    a[12] = (double) .0759;
    a[13] = (double) .0651;
    a[14] = (double) .0546;
    a[15] = (double) .0444;
    a[16] = (double) .0343;
    a[17] = (double) .0244;
    a[18] = (double) .0146;
    a[19] = (double) .0049;
  }
  else if (n == 41)
  {
    a[0] = (double) .394;
    a[1] = (double) .2719;
    a[2] = (double) .2357;
    a[3] = (double) .2091;
    a[4] = (double) .1876;
    a[5] = (double) .1693;
    a[6] = (double) .1531;
    a[7] = (double) .1384;
    a[8] = (double) .1249;
    a[9] = (double) .1123;
    a[10] = (double) .1004;
    a[11] = (double) .0891;
    a[12] = (double) .0782;
    a[13] = (double) .0677;
    a[14] = (double) .0575;
    a[15] = (double) .0476;
    a[16] = (double) .0379;
    a[17] = (double) .0283;
    a[18] = (double) .0188;
    a[19] = (double) .0094;
  }
  else if (n == 42)
  {
    a[0] = (double) .3917;
    a[1] = (double) .2701;
    a[2] = (double) .2345;
    a[3] = (double) .2085;
    a[4] = (double) .1874;
    a[5] = (double) .1694;
    a[6] = (double) .1535;
    a[7] = (double) .1392;
    a[8] = (double) .1259;
    a[9] = (double) .1136;
    a[10] = (double) .102;
    a[11] = (double) .0909;
    a[12] = (double) .0804;
    a[13] = (double) .0701;
    a[14] = (double) .0602;
    a[15] = (double) .0506;
    a[16] = (double) .0411;
    a[17] = (double) .0318;
    a[18] = (double) .0227;
    a[19] = (double) .0136;
    a[20] = (double) .0045;
  }
  else if (n == 43)
  {
    a[0] = (double) .3894;
    a[1] = (double) .2684;
    a[2] = (double) .2334;
    a[3] = (double) .2078;
    a[4] = (double) .1871;
    a[5] = (double) .1695;
    a[6] = (double) .1539;
    a[7] = (double) .1398;
    a[8] = (double) .1269;
    a[9] = (double) .1149;
    a[10] = (double) .1035;
    a[11] = (double) .0927;
    a[12] = (double) .0824;
    a[13] = (double) .0724;
    a[14] = (double) .0628;
    a[15] = (double) .0534;
    a[16] = (double) .0442;
    a[17] = (double) .0352;
    a[18] = (double) .0263;
    a[19] = (double) .0175;
    a[20] = (double) .0087;
  }
  else if (n == 44)
  {
    a[0] = (double) .3872;
    a[1] = (double) .2667;
    a[2] = (double) .2323;
    a[3] = (double) .2072;
    a[4] = (double) .1868;
    a[5] = (double) .1695;
    a[6] = (double) .1542;
    a[7] = (double) .1405;
    a[8] = (double) .1278;
    a[9] = (double) .116;
    a[10] = (double) .1049;
    a[11] = (double) .0943;
    a[12] = (double) .0842;
    a[13] = (double) .0745;
    a[14] = (double) .0651;
    a[15] = (double) .056;
    a[16] = (double) .0471;
    a[17] = (double) .0383;
    a[18] = (double) .0296;
    a[19] = (double) .0211;
    a[20] = (double) .0126;
    a[21] = (double) .0042;
  }
  else if (n == 45)
  {
    a[0] = (double) .385;
    a[1] = (double) .2651;
    a[2] = (double) .2313;
    a[3] = (double) .2065;
    a[4] = (double) .1865;
    a[5] = (double) .1695;
    a[6] = (double) .1545;
    a[7] = (double) .141;
    a[8] = (double) .1286;
    a[9] = (double) .117;
    a[10] = (double) .1062;
    a[11] = (double) .0959;
    a[12] = (double) .086;
    a[13] = (double) .0765;
    a[14] = (double) .0673;
    a[15] = (double) .0584;
    a[16] = (double) .0497;
    a[17] = (double) .0412;
    a[18] = (double) .0328;
    a[19] = (double) .0245;
    a[20] = (double) .0163;
    a[21] = (double) .0081;
  }
  else if (n == 46)
  {
    a[0] = (double) .383;
    a[1] = (double) .2635;
    a[2] = (double) .2302;
    a[3] = (double) .2058;
    a[4] = (double) .1862;
    a[5] = (double) .1695;
    a[6] = (double) .1548;
    a[7] = (double) .1415;
    a[8] = (double) .1293;
    a[9] = (double) .118;
    a[10] = (double) .1073;
    a[11] = (double) .0972;
    a[12] = (double) .0876;
    a[13] = (double) .0783;
    a[14] = (double) .0694;
    a[15] = (double) .0607;
    a[16] = (double) .0522;
    a[17] = (double) .0439;
    a[18] = (double) .0357;
    a[19] = (double) .0277;
    a[20] = (double) .0197;
    a[21] = (double) .0118;
    a[22] = (double) .0039;
  }
  else if (n == 47)
  {
    a[0] = (double) .3808;
    a[1] = (double) .262;
    a[2] = (double) .2291;
    a[3] = (double) .2052;
    a[4] = (double) .1859;
    a[5] = (double) .1695;
    a[6] = (double) .155;
    a[7] = (double) .142;
    a[8] = (double) .13;
    a[9] = (double) .1189;
    a[10] = (double) .1085;
    a[11] = (double) .0986;
    a[12] = (double) .0892;
    a[13] = (double) .0801;
    a[14] = (double) .0713;
    a[15] = (double) .0628;
    a[16] = (double) .0546;
    a[17] = (double) .0465;
    a[18] = (double) .0385;
    a[19] = (double) .0307;
    a[20] = (double) .0229;
    a[21] = (double) .0153;
    a[22] = (double) .0076;
  }
  else if (n == 48)
  {
    a[0] = (double) .3789;
    a[1] = (double) .2604;
    a[2] = (double) .2281;
    a[3] = (double) .2045;
    a[4] = (double) .1855;
    a[5] = (double) .1693;
    a[6] = (double) .1551;
    a[7] = (double) .1423;
    a[8] = (double) .1306;
    a[9] = (double) .1197;
    a[10] = (double) .1095;
    a[11] = (double) .0998;
    a[12] = (double) .0906;
    a[13] = (double) .0817;
    a[14] = (double) .0731;
    a[15] = (double) .0648;
    a[16] = (double) .0568;
    a[17] = (double) .0489;
    a[18] = (double) .0411;
    a[19] = (double) .0335;
    a[20] = (double) .0259;
    a[21] = (double) .0185;
    a[22] = (double) .0111;
    a[23] = (double) .0037;
  }
  else if (n == 49)
  {
    a[0] = (double) .377;
    a[1] = (double) .2589;
    a[2] = (double) .2271;
    a[3] = (double) .2038;
    a[4] = (double) .1851;
    a[5] = (double) .1692;
    a[6] = (double) .1553;
    a[7] = (double) .1427;
    a[8] = (double) .1312;
    a[9] = (double) .1205;
    a[10] = (double) .1105;
    a[11] = (double) .101;
    a[12] = (double) .0919;
    a[13] = (double) .0832;
    a[14] = (double) .0748;
    a[15] = (double) .0667;
    a[16] = (double) .0588;
    a[17] = (double) .0511;
    a[18] = (double) .0436;
    a[19] = (double) .0361;
    a[20] = (double) .0288;
    a[21] = (double) .0215;
    a[22] = (double) .0143;
    a[23] = (double) .0071;
  }
  else if (n == 50)
  {
    a[0] = (double) .3751;
    a[1] = (double) .2574;
    a[2] = (double) .226;
    a[3] = (double) .2032;
    a[4] = (double) .1847;
    a[5] = (double) .1691;
    a[6] = (double) .1554;
    a[7] = (double) .143;
    a[8] = (double) .1317;
    a[9] = (double) .1212;
    a[10] = (double) .1113;
    a[11] = (double) .102;
    a[12] = (double) .0932;
    a[13] = (double) .0846;
    a[14] = (double) .0764;
    a[15] = (double) .0685;
    a[16] = (double) .0608;
    a[17] = (double) .0532;
    a[18] = (double) .0459;
    a[19] = (double) .0386;
    a[20] = (double) .0314;
    a[21] = (double) .0244;
    a[22] = (double) .0174;
    a[23] = (double) .0104;
    a[24] = (double) .0035;
  }

  if (n > 50 || n < 3)
  {
#ifdef NOISY
    printf ("  THIS IS THE SHAPIRO-WILK TEST FOR SMALL SAMPLES\n");
    printf ("  THE SAMPLE SIZE MUST BE LESS THAN OR EQUAL TO 50\n");
#endif				/* NOISY */
    y[0] = y[1] = 0.0;
  }
  else
  {
    /* following bug reported by Anders Ledberg <anders@pet.neuro.ki.se>
       on 8 May 1996 */
    /*-for (i = 1; i <= k; ++i) 
       sumb += a[i - 1] * (x[n - i + 1] - x[i]);  */
    for (i = 0; i < k; ++i)
      sumb += a[i] * (xcopy[n - i + 1] - xcopy[i]);
    y[0] = sumb * sumb / s2;
    y[1] = s2;

#ifdef NOISY
    printf ("  TEST13 SW(N)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  }
  free (xcopy);
  return y;
}

/* this is actually the Weisberg-Bingham stat. I need to
OCR the constants and implment this correctly */

double *shapiro_francia ( double * x,int n)

{
  static double y[2];
  double suma = 0.0, sumb = 0.0, sumc = 0.0, sumd = 0.0, z, *xcopy;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in shapiro_francia\n"), exit (-1);

  for (i = 0; i < n; ++i)
    xcopy[i] = x[i];
  qsort (xcopy, n, sizeof (double), dcmp);

  for (i = 0; i < n; ++i)
  {
    z = xinormal ((i + 1 - 0.375) / (n + 0.25));
    suma += z * xcopy[i];	
    sumb += z * z;
    sumc += xcopy[i];
    sumd += xcopy[i] * xcopy[i];
  }
  y[0] = suma * suma / sumb / (sumd - sumc * sumc / n);
#ifdef NOISY
  printf ("  TEST14 SF(N)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  free (xcopy);
  return y;
}				/* test14_ */
#include<stdio.h>
#include<math.h>

double *shapiro_wilk_exp (double * x,int n)

{
  static double y[2];
  double mean, b, s1, xs, sum1 = 0.0, sum2 = 0.0;
  int i;

  for (i = 0; i < n; ++i)
    if (i == 0 || xs > x[i])
      xs = x[i];

  for (i = 0; i < n; ++i)
  {
    sum1 += x[i];
    sum2 += x[i] * x[i];
  }
  s1 = sum2 - sum1 * sum1 / n;
  mean=sum1/n;
  b = (mean - xs) * sqrt ((double) n / (n - 1.0));
  y[0] = b * b / s1;

#ifdef NOISY
  printf ("  TEST15 SW(E)  =%10.4f\n", y[0]);
#endif				/* NOISY */

  return y;
}

double *watson_u2_exp (double *x,   int n)
{
  double *xcopy, mean = 0.0, zbar = 0.0, fn2, fx, sum4 = 0.0;
  static double y[2];
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in watson_u2_exp\n"), exit (-1);

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
  }
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);
  for (i = 0; i < n; ++i)
  {
    fx = 1 - exp (-xcopy[i] / mean);
    if (fx <= 1e-5)
      fx = 1e-5;
    if (fx >= 0.99999)
      fx = 0.99999;
    /* sum3 += 2 * (i + 1) * (log (fx) + log (1.0 - fx[n - i - 1])); */
    fn2 = (2.0 * i + 1.0) / (2.0 * n);
      sum4 += (fx - fn2) * (fx - fn2);
    fn2 = (2.0 * (i + 1) - 1.0) / (2.0 * n);
    zbar += fx;
  }
  zbar /= n;
  y[0] = (1.0 / (n * 12) + sum4) - n * (zbar - .5) * (zbar - .5);
  y[0] *= 1.0 + 0.16 / n;
#ifdef NOISY
  printf ("  TEST19 WU2(E) =%10.4f\n", y[0]);
#endif				/* NOISY */
  free (xcopy);
  return y;
}

double *watson_u2 (double * x,int  n)

{
  double *xcopy, mean = 0.0, sdx = 0.0, sqrt2, zbar = 0.0; 
  double fn2, fx, sum4 = 0.0;
  static double y[2];
  int i;

  sqrt2 = sqrt ((double) 2.0);

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in anderson_darling\n"), exit (-1);
 
  for (i = 0; i < n; ++i)
  {
    xcopy[i] = x[i];
    mean += x[i];
    sdx += x[i] * x[i];
  }
  sdx = sqrt ((n * sdx - mean * mean) / (n * (n - 1)));
  mean /= n;
  qsort (xcopy, n, sizeof (double), dcmp);

  for (i = 0; i < n; ++i)
  {
    xcopy[i] = (xcopy[i] - mean) / sdx;
    fn2 = (2.0 * (i + 1) - 1.0) / (2.0 * n);
    fx= 0.5 + normp (xcopy[i] / sqrt2) / 2.0;
    if (fx <= 0.0)
      fx = 1e-5;
    if (fx >= 1.0)
      fx = 0.99999;
    zbar += fx;
    sum4 += (fx - fn2) * (fx - fn2);
  }

  zbar /= n;
  y[0] = (1.0 / (n * 12) + sum4) - n * (zbar - .5) * (zbar - .5);
  y[0] *= 0.5 / n + 1.0;

#ifdef NOISY
  printf ("  TEST6  WU2(N) =%10.4f\n", y[0]);
#endif				/* NOISY */
  free(xcopy);
  return y;
}

double *weisberg_bingham (double * x,int n)

{
  static double y[2];
  double suma = 0.0, sumb = 0.0, sumc = 0.0, sumd = 0.0, z, *xcopy;
  int i;

  if ((xcopy = (double *) malloc (n * sizeof (double))) == NULL)
    fprintf (stderr, "Memory error in shapiro_francia\n"), exit (-1);

  for (i = 0; i < n; ++i)
    xcopy[i] = x[i];
  qsort (xcopy, n, sizeof (double), dcmp);

  for (i = 0; i < n; ++i)
  {
    z = xinormal ((i + 1 - 0.375) / (n + 0.25));
    suma += z * xcopy[i];	
    sumb += z * z;
    sumc += xcopy[i];
    sumd += xcopy[i] * xcopy[i];
  }
  y[0] = suma * suma / sumb / (sumd - sumc * sumc / n);
#ifdef NOISY
  printf ("  TEST14 SF(N)  =%10.4f\n", y[0]);
#endif				/* NOISY */
  free (xcopy);
  return y;
}				/* test14_ */

double xinormal (double pee)

{
  double f0, pind, pw, px;
  static double p[5] = {-0.322232431088, -1., -0.342242088547,
  -0.0204231210245, -4.53642210148e-5};
  static double q[5] = {0.099348462606, 0.588581570495, 0.531103462366,
  0.10353775285, 0.0038560700634};

  pind = pee;
  if (pee < 1e-10)
    return (double) -10.0;
  else if (pee >= 1.0)
    return (double) 10.0;
  else if (pee == 0.5)
    return (double) 0.5;
  else
  {
    if (pee > .5)
      pee--;
    pw = sqrt (log (1 / (pee * pee)));
    f0 = (((pw * q[4] + q[3]) * pw + q[2]) * pw + q[1]) * pw + q[0];
    px = pw + ((((pw * p[4] + p[3]) * pw + p[2])
		* pw + p[1]) * pw + p[0]) / f0;
    if (pind < .5)
      px = -px;
    return px;
  }
}

double * coeff_variation(double * x,int n)

{
    static double y[2];
    double s2, s3, s4, sum2=0.0, sum4=0.0;
    int i;

    for (i = 0; i < n; ++i) 
	sum2 += log(x[i]);
    for (i = 0; i < n; ++i) 
	sum4 += (log(x[i]) - sum2 / n) * (log(x[i]) - sum2 / n);
    s2 = sum4 / (n - 1);
    s3 = exp(s2) - 1;
    s4 = sqrt(s3);
    y[0] = sqrt(exp(s4/(n-1))-1);
#ifdef NOISY
 printf("  TEST23 CV(L)  =%10.4f\n",s4);
#endif /* NOISY */
    return y;
}

double *mod_maxlik_ratio (double * x,int  n)

{
  static double y[2];
  double mean = 0.0;
  double e1, m1, s1, s2, s3;
  double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0;
  int i;

  for (i = 0; i < n; ++i)
  {
    sum1 += x[i];
    sum2 += log (x[i]);
  }
  s1 = (n * sum2 - mean * mean) / (n * (n - 1));
  mean = sum1 / n;
  m1 = sum2 / n;

  for (i = 0; i < n; ++i)
  {
    sum4 += (log (x[i]) - m1) * (log (x[i]) - m1);
    sum3 += pow (x[i] - mean, 3.0);
  }
  s2 = sum4 / n;

  if (sum3 < 0.0)
  {
#ifdef NOISY
    for (i = 0; i < n; ++i)
    {
      printf (" %g", x[i]);
      if (i % 4)
	printf ("\n");
    }
    printf (" THIRD SAMPLE MOMENT ABOUT THE MEAN IS LESS THAN ZERO\n");
    printf (" HENCE WE ACCEPT THE NULL HYPOTHESIS OF NORMALITY\n");
#endif /* NOISY */
    y[0] = 0.0;
  }
  else
  {
    for (i = 0; i < n; ++i)
      sum5 += (x[i] - mean) * (x[i] - mean);
    s1 = sqrt (sum5 / n);
    s3 = sqrt (s2);
    e1 = exp (m1);
    y[0] = s1 / (s3 * e1);
  }
#ifdef NOISY
  printf ("  TEST22 LR(NL) =%10.4f\n", y[0]);
#endif				/* NOISY */
  return y;
}

