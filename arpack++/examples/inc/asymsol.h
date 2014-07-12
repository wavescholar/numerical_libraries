/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ASymSol.cc
   Template functions that exemplify how to print information
   about symmetric standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ASYMSOL_H
#define ASYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "arlsmat.h"

template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnz, FLOAT A[], INT irow[], INT pcol[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/

{

  INT                  i;
  FLOAT*               Ax;
  FLOAT*               ResNorm;
  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x \n \n";

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.


template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnzA, FLOAT A[], INT irowA[],
              INT pcolA[], INT nnzB, FLOAT B[], INT irowB[], INT pcolB[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of symmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  INT                  i;
  FLOAT                *Ax, *Bx;
  FLOAT                *ResNorm;
  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl << endl;

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    Bx      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution.


#endif // ASYMSOL_H

