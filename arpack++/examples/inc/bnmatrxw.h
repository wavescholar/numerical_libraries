/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BNMatrxW.h
   Class template for the n by n matrix A with entries

          A(i,j) = k*(s(i))*(t(j)-1), if j-nl<= i <= j;
                   k*(t(j))*(s(i)-1), if j < i <= j+nu,

   where s(i) = i/(n+1), t(j) = j/(n+1) and k = 1/(n+1).       

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef BNMATRXW_H
#define BNMATRXW_H


template<class FLOAT>
void MatrixW(int n, int nl, int nu, FLOAT* &A)
{

  int   i, j, p;
  FLOAT h;

  // Defining h.

  h = 1.0/FLOAT(n+1);

  // Creating output vector A.

  A = new FLOAT[n*(nl+nu+1)];

  p = 0;

  for (j=1; j<=n; j++) {

    for (i=j-nu; i<=j; i++) {
      if (i<1) {
        A[p++] = 0.0;
      }
      else {
        A[p++] = h*(FLOAT(j)*h-1.0)*(FLOAT(i)*h);
      }
    }

    for(i=j+1; i<=(j+nl); i++) {
      if (i>n) {
        A[p++] = 0.0;
      }
      else {
        A[p++]  = h*(FLOAT(i)*h-1.0)*(FLOAT(j)*h);
      }
    }

  }

} // MatrixW.

#endif // BNMATRXW_H

