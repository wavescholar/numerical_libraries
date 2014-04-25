/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BCMatrxB.h
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

#ifndef BCMATRXB_H
#define BCMATRXB_H

#include "arcomp.h"

template<class FLOAT, class INT>
void CompMatrixB(INT n, arcomplex<FLOAT> rho, INT& nL, 
                 INT& nU, arcomplex<FLOAT>* &A)
{

  INT              i;
  arcomplex<FLOAT> dd, dl, du, h;

  // Defining constants.

  const arcomplex<FLOAT> one( 1.0, 0.0);
  const arcomplex<FLOAT> two( 2.0, 0.0);
  const arcomplex<FLOAT> half( 0.5, 0.0);

  h  = one/arcomplex<FLOAT>((FLOAT)(n+1),0.0);
  dd = two/h;
  dl = -(one/h) - half*rho;
  du = -(one/h) + half*rho;

  // Defining the lower and the upper bandwidth.

  nL = 1;
  nU = 1;

  // Creating output vector A.

  A  = new arcomplex<FLOAT>[3*n];
  for (i=1; i<(3*n); i+=3) {
    if (i-1)   A[i-1] = du;
               A[i]   = dd;
    if (n-i-1) A[i+1] = dl;
  }

} //  CompMatrixB.


#endif // BCMATRXB_H

