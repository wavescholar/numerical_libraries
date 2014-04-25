/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LCMatrxB.h
   Function template for the tridiagonal matrix derived from 
   the standard central difference of the 1-d convection diffusion 
   operator u" + rho*u' on the interval [0, 1] with zero
   Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXB_H
#define LCMATRXB_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixB(INT n, arcomplex<FLOAT> rho, INT& nnz, 
                 arcomplex<FLOAT>* &A, INT* &irow, INT* &pcol) 
{

  INT              i, j;
  arcomplex<FLOAT> dd, dl, du, s, h, h2;

  // Defining constants.

  const arcomplex<FLOAT> one( 1.0, 0.0);
  const arcomplex<FLOAT> two( 2.0, 0.0);

  h  = one/arcomplex<FLOAT>((FLOAT)(n+1),0.0);
  h2 = h*h;
  s  = rho/two;
  dd = two/h2;
  dl = -(one/h2) - s/h;
  du = -(one/h2) + s/h;

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
      A[j++]  = du;
    }

    irow[j]   = i;
    A[j++]    = dd;

    if (i!=(n-1)) {
      irow[j] = i+1;
      A[j++]  = dl;
    }

    pcol[i+1] = j;

  }

} //  CompMatrixB.


#endif // LCMATRXB_H

