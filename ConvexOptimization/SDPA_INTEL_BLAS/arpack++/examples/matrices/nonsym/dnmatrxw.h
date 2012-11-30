/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE DNMatrxW.h
   Class template for the n by n matrix A with entries

          A(i,j) = k*(s(i))*(t(j)-1), if i <= j,
                   k*(t(j))*(s(i)-1), if j < i.

   where s(i) = i/(n+1), t(j) = j/(n+1) and k = 1/(n+1).       

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXW_H
#define DNMATRXW_H


template<class FLOAT>
void MatrixW(int m, int n, FLOAT* &A)
{

  int   i, j, p;
  FLOAT h, k;

  // Defining h.

  h = 1.0/FLOAT(m+1);
  k = 1.0/FLOAT(n+1);

  // Creating output vector A.

  A = new FLOAT[m*n];

  p = 0;
  for (j=1; j<=n; j++) {

    for (i=1; i<=j; i++) {
      A[p++] = k*(FLOAT(j)*k-1.0)*(FLOAT(i)*h);
    }
    for(i=j+1; i<=m; i++) {
      A[p++]  = k*(FLOAT(i)*h-1.0)*(FLOAT(j)*k);
    }

  }

} // MatrixW.

#endif // DNMATRXW_H

