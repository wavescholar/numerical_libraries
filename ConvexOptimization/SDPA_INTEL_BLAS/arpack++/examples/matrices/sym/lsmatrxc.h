/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LSMatrxC.h
   Function template for the one dimensional discrete Laplacian
   on the interval [0, 1], with zero Dirichlet boundary conditions.
   The difference between the function shown here and SymmetricMatrixB
   is only the scaling factor.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LSMATRXC_H
#define LSMATRXC_H

#include <math.h>

template<class FLOAT, class INT>
void SymmetricMatrixC(INT n, INT& nnz, FLOAT* &A, 
                      INT* &irow, INT* &pcol, char uplo = 'L')

{

  // Defining internal variables.

  INT    i, j;
  FLOAT  h, df, dd;

  // Defining constants.

  h  = FLOAT(n+1);
  dd = 2.0*h;
  df = -h;

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
} // SymmetricMatrixC.

#endif // LSMATRXC_H

