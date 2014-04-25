/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE DCMatrxB.h
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

#ifndef DCMATRXB_H
#define DCMATRXB_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixB(INT n, arcomplex<FLOAT> rho, arcomplex<FLOAT>* &A)
{

  INT              i, j;
  arcomplex<FLOAT> dd, dl, du, h;

  // Defining constants.

  const arcomplex<FLOAT> one( 1.0, 0.0);
  const arcomplex<FLOAT> two( 2.0, 0.0);
  const arcomplex<FLOAT> half( 0.5, 0.0);

  h  = one/arcomplex<FLOAT>((FLOAT)(n+1),0.0);
  dd = two/h;
  dl = -(one/h) - half*rho;
  du = -(one/h) + half*rho;

  // Creating output vector A.

  A  = new arcomplex<FLOAT>[n*n];

  for (i=0; i<n*n; i++) A[i] = arcomplex<FLOAT>(0.0, 0.0);

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)      A[j-1]  = du;
                A[j]    = dd;
    if (n-i-1)  A[j+1]  = dl;
  }

} //  CompMatrixB.


#endif // DCMATRXB_H

