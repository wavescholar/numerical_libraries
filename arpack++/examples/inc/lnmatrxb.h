/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE LNMatrxB.h
   Function template for the matrix

                          | T -I          |
                          |-I  T -I       |
                      A = |   -I  T       |
                          |        ...  -I|
                          |           -I T|

   derived from the standard central difference discretization of the
   2-dimensional convection-diffusion operator (Laplacian u) + rho*(du/dx)
   on a unit square with zero Dirichlet boundary conditions.
   When rho*h/2 <= 1, the discrete convection-diffusion operator has real
   eigenvalues.  When rho*h/2 > 1, it has COMPLEX eigenvalues.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LNMATRXB_H
#define LNMATRXB_H

#include <math.h>

template<class FLOAT, class INT>
void BlockTridMatrix(INT nx, INT& n, INT& nnz,
                     FLOAT* &A, INT* &irow, INT* &pcol)

{

  // Defining internal variables.

  INT    i, j;
  FLOAT  h, h2, df;
  FLOAT  dd, dl, du;

  // Defining constants.

  const FLOAT rho  = 0.0;

  h   = 1.0/FLOAT(nx+1);
  h2  = h*h;
  dd  = 4.0/h2;
  df  = -1.0/h2;
  dl  = df - 5.0e-1*rho/h;
  du  = df + 5.0e-1*rho/h;

  // Defining the number of columns and nonzero elements of matrix.

  n   = nx*nx;
  nnz = 5*n-4*nx;

  // Creating output vectors.

  A    = new FLOAT[nnz];
  irow = new INT[nnz];
  pcol = new INT[n+1];

  // Creating matrix.

  pcol[0] = 0;
  i       = 0;

  for (j = 0; j < n; j++) {
    if (j >= nx) {
      A[i] = df;   irow[i++] = j-nx;
    }
    if ((j%nx) != 0) {
      A[i] = du;   irow[i++] = j-1;
    }
    A[i] = dd;     irow[i++] = j;
    if (((j+1)%nx) != 0) {
      A[i] = dl;   irow[i++] = j+1;
    }
    if (j < n-nx) {
      A[i] = df;   irow[i++] = j+nx;
    }
    pcol[j+1] = i;
  }

} // BlockTridMatrix.

#endif // LNMATRXB_H

