/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE DCMatrxC.h
   Function template for the mass matrix formed by using piecewise 
   linear elements on [0, 1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DCMATRXC_H
#define DCMATRXC_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixC(INT n, arcomplex<FLOAT>* &A)
{

  INT              i, j;
  arcomplex<FLOAT> diag, sub;

  // Defining constants.

  sub  = arcomplex<FLOAT>(1.0/(FLOAT)(n+1), 0.0);
  diag = arcomplex<FLOAT>(4.0/(FLOAT)(n+1), 0.0);

  // Creating output vector A.

  A  = new arcomplex<FLOAT>[n*n];

  for (i=0; i<n*n; i++) A[i] = arcomplex<FLOAT>(0.0, 0.0);

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)      A[j-1]  = sub;
                A[j]    = diag;
    if (n-i-1)  A[j+1]  = sub;
  }

} // CompMatrixC.

#endif // DCMATRXC_H

