/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LCMatrxA.h
   Function template for the nx*nx by nx*nx block tridiagonal matrix

                              | T -I          |
                              |-I  T -I       |
                         OP = |   -I  T       |
                              |        ...  -I|
                              |           -I T|

   derived from the standard central difference discretization
   of the 2 dimensional convection-diffusion operator
                       (Laplacian u) + rho*(du/dx) 
   on a unit square with zero boundary conditions.
   T is a nx by nx tridiagonal matrix with DD on the diagonal,
   DL on the subdiagonal, and DU on the superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LCMATRXA_H
#define LCMATRXA_H

#include "arcomp.h"
#include "blas1c.h"

template<class FLOAT, class INT>
void CompMatrixA(INT nx, INT& nnz, arcomplex<FLOAT>* &A, 
                 INT* &irow, INT* &pcol)
{

  INT              i, j, k, id;
  arcomplex<FLOAT> h, h2, dd, dl, du, f;

  // Defining constants.

  const arcomplex<FLOAT> half(0.5,0.0);
  const arcomplex<FLOAT> one(1.0,0.0);
  const arcomplex<FLOAT> four(4.0,0.0);
  const arcomplex<FLOAT> rho(1.0e2,0.0);

  h   = one/arcomplex<FLOAT>((FLOAT)(nx+1),0);
  h2  = h*h;
  f   = -(one/h2);
  dd  = four/h2;
  dl  = f - half*rho/h;
  du  = f + half*rho/h;

  // Defining the number of nonzero matrix elements.

  nnz = (5*nx-4)*nx;

  // Creating output vectors.

  A    = new arcomplex<FLOAT>[nnz];
  irow = new INT[nnz];
  pcol = new INT[nx*nx+1];

  // Filling A, irow and pcol.

  pcol[0] = 0;
  j       = 0;
  id      = 0;

  for (k=0; k!=nx; k++) {
    for (i=0; i!=nx; i++) {

      if (k) {
        irow[j] = id-nx;
        A[j++]  = f; 
      }

      if (i) {
        irow[j] = id-1;
        A[j++]  = du;
      }

      irow[j]   = id;
      A[j++]    = dd;

      if (i!=(nx-1)) {
        irow[j] = id+1;
        A[j++]  = dl;
      }

      if (k!=(nx-1)) {
        irow[j] = id+nx;
        A[j++]  = f;
      }

      pcol[++id]= j;
    }
  }     
  
} // CompMatrixA.


#endif // LCMATRXA_H
