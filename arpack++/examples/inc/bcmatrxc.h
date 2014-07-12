/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BCMatrxC.h
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

#ifndef BCMATRXC_H
#define BCMATRXC_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixC(INT n, INT& nL, INT& nU, arcomplex<FLOAT>* &A)
{

  INT              i;
  arcomplex<FLOAT> diag, sub;

  // Defining constants.

  sub  = arcomplex<FLOAT>(1.0/(FLOAT)(n+1), 0.0);
  diag = arcomplex<FLOAT>(4.0/(FLOAT)(n+1), 0.0);

  // Defining the lower and upper bandwidth.

  nL = 1;
  nU = 1;

  // Creating output vector A.

  A  = new arcomplex<FLOAT>[3*n];
  for (i=1; i<(3*n); i+=3) {
    if (i-1)   A[i-1] = sub;
               A[i]   = diag;
    if (n-i-1) A[i+1] = sub;
  }

} // CompMatrixC.

#endif // BCMATRXC_H

