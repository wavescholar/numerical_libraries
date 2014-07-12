/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LNMatrxF.h
   Funtion template that generates a tridiagonal nonsymmetric matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXF_H
#define LNMATRXF_H


template<class FLOAT, class INT>
void NonSymMatrixF(INT n, INT& nnz, FLOAT* &A, INT* &irow, INT* &pcol)
{

  INT   i,j;
  FLOAT diag, sub;

  // Defining constants.

  sub  = 1.0;
  diag = 4.0;

  // Defining the number of nonzero matrix elements.

  nnz = 3*n-2;

  // Creating output vectors.

  A    = new FLOAT[nnz];
  irow = new INT[nnz];
  pcol = new INT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j = 0;
  for (i=0; i!=n; i++) {
    if (i != 0) {
      irow[j] = i-1;
      A[j++]  = sub;
    }
    irow[j] = i;
    A[j++]  = diag;
    if (i != (n-1)) {
      irow[j] = i+1;
      A[j++]  = sub;
    }
    pcol[i+1] = j;
  }

} // NonSymMatrixF.


#endif // LNMATRXF_H

