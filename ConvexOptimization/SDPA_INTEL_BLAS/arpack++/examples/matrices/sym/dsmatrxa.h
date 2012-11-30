/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE DSMatrxA.h
   Function template for the matrix

                          | T -I          |
                          |-I  T -I       |
                      A = |   -I  T       |
                          |        ...  -I|
                          |           -I T|

   derived from the standard central difference discretization of the
   2-dimensional Laplacian on the unit square with zero Dirichlet 
   boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef DSMATRXA_H
#define DSMATRXA_H

#include <math.h>

template<class FLOAT, class INT>
void DenseMatrixA(INT nx, INT& n, FLOAT* &A, char uplo = 'L')
{

  // Defining internal variables.

  INT    i, j;
  FLOAT  h2, df, dd;

  // Defining constants.

  h2  = 1.0/(FLOAT(nx+1)*FLOAT(nx+1));
  dd  = 4.0/h2;
  df  = -1.0/h2;

  // Defining the dimension of A.

  n   = nx*nx;

  // Creating matrix A.

  A   = new FLOAT[(n*n+n)/2];

  for (i=0; i<(n*n+n)/2; i++) A[i]=0.0;

  if (uplo == 'L') {

    for (i=0, j=0; i<n; j+=(n-(i++))) {

      A[j] = dd;
      if (((i+1)%nx) != 0) {
        A[j+1] = df;
      }
      if (i < n-nx) {
        A[j+nx] = df;
      }
    
    }

  }
  else { // uplo = "U"

    for (i=0, j=0; i<n; j+=(++i)) {

      A[j+i] = dd;
      if ((i%nx) != 0) {
        A[j+i-1] = df;
      }
      if (i >= nx) {
        A[j+i-nx] = df;
      }
    
    }  

  }

} // DenseMatrixA.

#endif // DSMATRXA_H

