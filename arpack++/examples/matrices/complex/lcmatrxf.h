/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LCMatrxF.h
   Function template for a tridiagonal complex matrix.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXF_H
#define LCMATRXF_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixF(INT n, INT& nnz, arcomplex<FLOAT>* &A, 
                 INT* &irow, INT* &pcol)
{

  INT              i, j;
  arcomplex<FLOAT> h, dd, ds;

  // Defining constants.

  const arcomplex<FLOAT> one(1.0, 0.0);
  const arcomplex<FLOAT> four(4.0, 0.0);

  h  = one/arcomplex<FLOAT>((FLOAT)(n+1),0.0);
  dd = four*h;
  ds = one*h;

  // Defining the number of nonzero matrix elements.

  nnz = 3*n-2;

  // Creating output vectors.

  A    = new arcomplex<FLOAT>[nnz];
  irow = new INT[nnz];
  pcol = new INT[n+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j       = 0;

  for (i=0; i!=n; i++) {

    if (i) {
      irow[j] = i-1;
      A[j++]  = ds;
    }

    irow[j]   = i;
    A[j++]    = dd;

    if (i!=(n-1)) {
      irow[j] = i+1;
      A[j++]  = ds;
    }

    pcol[i+1] = j;

  }

} // CompMatrixF.


#endif // LCMATRXF_H

