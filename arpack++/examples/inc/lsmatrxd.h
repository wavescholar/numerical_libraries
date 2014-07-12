/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LSMatrxD.h
   Function template for the 1-dimensional mass matrix
   on the interval [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LSMATRXD_H
#define LSMATRXD_H

#include <math.h>

template<class FLOAT, class INT>
void SymmetricMatrixD(INT n, INT& nnz, FLOAT* &A, 
                      INT* &irow, INT* &pcol, char uplo = 'L')

{

  // Defining internal variables.

  INT    i, j;
  FLOAT  h, df, dd;

  // Defining constants.

  h  = 1.0/FLOAT(n+1);
  dd = 4.0/6.0*h;
  df = h/6.0;

  // Defining the number of nonzero elements in A.

  nnz = 2*n-1;

  // Creating output vectors.

  A    = new FLOAT[nnz];
  irow = new INT[nnz];
  pcol = new INT[n+1];

  // Defining matrix A.

  pcol[0] = 0;
  i       = 0;

  if (uplo == 'U') {

    for (j = 0; j < n; j++) {
      if (j) {
        A[i] = df;   irow[i++] = j-1;
      }
      A[i] = dd;     irow[i++] = j;
      pcol[j+1] = i;
    }

  }
  else {

    for (j = 0; j < n; j++) {
      A[i] = dd;     irow[i++] = j;
      if (n-j-1) {
        A[i] = df;   irow[i++] = j+1;
      }
      pcol[j+1] = i;
    }

  }
} // SymmetricMatrixD.

#endif // LSMATRXD_H

