/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LNMatrxE.h
   Function template that generates a nonsymmetric 
   tridiagonal matrix with 2 on the main diagonal,
   3 on the superdiagonal and -2 on the subdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXE_H
#define LNMATRXE_H


template<class FLOAT, class INT>
void NonSymMatrixE(INT n, INT& nnz, FLOAT* &A, INT* &irow, INT* &pcol)
{

  INT   i, j;

  // Defining constants.
  
  const FLOAT three = 3.0;
  const FLOAT two   = 2.0;

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
      A[j++]  = three;
    }
    irow[j] = i;
    A[j++]  = two;
    if (i != (n-1)) {
      irow[j] = i+1;
      A[j++]  = -two;
    }
    pcol[i+1] = j;
  }

} // NonSymMatrixE.


#endif // LNMATRXE_H



