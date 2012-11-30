/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE DNMatrxB.h
   Function template for the stiffness matrix obtained from
   the finite element discretization of the 1-dimensional
   convection diffusion operator d^2u/dx^2 + rho*(du/dx) on
   the interval [0,1] with zero Dirichlet boundary conditions
   using linear elements.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DNMATRXB_H
#define DNMATRXB_H

template<class FLOAT, class INT>
void StiffnessMatrix(INT n, FLOAT rho, FLOAT* &A)
{

  INT   i, j;
  FLOAT dd, dl, du, h;

  // Defining constants.

  h  = 1.0/FLOAT(n+1);
  dd = 2.0/h;
  dl = -1.0/h - 0.5*rho;
  du = -1.0/h + 0.5*rho;

  // Creating output vector A.

  A  = new FLOAT[n*n];

  for (i=0; i<n*n; i++) A[i]=0.0;

  for (i=0, j=0; i<n; i++, j+=n+1) {
    if (i)     A[j-1] = du;
               A[j]   = dd;
    if (n-i-1) A[j+1] = dl;
  }

} //  StiffnessMatrix.

#endif // DNMATRXB_H
