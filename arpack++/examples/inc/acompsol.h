/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ACompSol.cc
   Template functions that exemplify how to print information
   about complex standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ACOMPSOL_H
#define ACOMPSOL_H

#include <math.h>
#include "arcomp.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arlnsmat.h"

template<class TYPE, class INT>
void Solution(INT nconv, INT n, INT nnz, TYPE A[], INT irow[],
              INT pcol[], TYPE EigVal[], TYPE* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problems
  on standard "cout" stream.
*/

{

  INT                    i;
  TYPE*                  Ax;
  double*                ResNorm;
  ARluNonSymMatrix<TYPE> matrix(n, nnz, A, irow, pcol);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "complex standard eigenvalue problem: A*x - lambda*x \n \n";

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

    Ax      = new TYPE[n];
    ResNorm = new double[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
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


template<class TYPE, class INT>
void Solution(INT nconv, INT n, INT nnzA, TYPE A[], INT irowA[],
              INT pcolA[], INT nnzB, TYPE B[], INT irowB[], INT pcolB[],
              TYPE EigVal[], TYPE* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  INT                     i;
  TYPE                   *Ax;
  TYPE                   *Bx;
  double                 *ResNorm;
  ARluNonSymMatrix<TYPE> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<TYPE> matrixB(n, nnzB, B, irowB, pcolB);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Complex generalized eigenvalue problem: A*x - lambda*B*x \n \n";

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

    Ax      = new TYPE[n];
    Bx      = new TYPE[n];
    ResNorm = new double[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/lapy2(real(EigVal[i]),imag(EigVal[i]));
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


#endif // ACOMPSOL_H

