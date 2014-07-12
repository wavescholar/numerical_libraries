/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ANSymSol.cc
   Template functions that exemplify how to print information
   about nonsymmetric standard and generalized eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ANSYMSOL_H
#define ANSYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "lapackc.h"
#include "arlnsmat.h"

template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnz, FLOAT A[], INT irow[], INT pcol[],
              FLOAT EigValR[], FLOAT EigValI[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problem
  on standard "cout" stream.
*/

{

  INT                     i;
  FLOAT*                  Ax;
  FLOAT*                  ResNorm;
  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real nonsymmetric eigenvalue problem: A*x - lambda*x \n \n";

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
    if (EigValI[i]>=0.0) {
      cout << " + " << EigValI[i] << " I" << endl;
    }
    else {
      cout << " - " << fabs(EigValI[i]) << " I" << endl;
    }
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (EigValI[i]==0.0) { // Eigenvalue is real.

        matrix.MultMv(&EigVec[i*n], Ax);
        axpy(n, -EigValR[i], &EigVec[i*n], 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigValR[i]);

      }
      else {                 // Eigenvalue is complex.

        matrix.MultMv(&EigVec[i*n], Ax);
        axpy(n, -EigValR[i], &EigVec[i*n], 1, Ax, 1);
        axpy(n, EigValI[i], &EigVec[(i+1)*n], 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        matrix.MultMv(&EigVec[(i+1)*n], Ax);
        axpy(n, -EigValI[i], &EigVec[i*n], 1, Ax, 1);
        axpy(n, -EigValR[i], &EigVec[(i+1)*n], 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(EigValR[i], EigValI[i]);
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
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
              FLOAT EigValR[], FLOAT EigValI[], FLOAT* EigVec = 0)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric generalized
  eigen-problem on standard "cout" stream.
*/

{

  INT                     i;
  FLOAT                   *Ax;
  FLOAT                   *Bx, *Bx1;
  FLOAT                   *ResNorm;
  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real nonsymmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl << endl;

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
    if (EigValI[i]>=0.0) {
      cout << " + " << EigValI[i] << " I" << endl;
    }
    else {
      cout << " - " << fabs(EigValI[i]) << " I" << endl;
    }
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    Bx      = new FLOAT[n];
    Bx1     = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (EigValI[i]==0.0) { // Eigenvalue is real.

        matrixA.MultMv(&EigVec[i*n], Ax);
        matrixB.MultMv(&EigVec[i*n], Bx);
        axpy(n, -EigValR[i], Bx, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigValR[i]);

      }
      else {                 // Eigenvalue is complex.

        matrixA.MultMv(&EigVec[i*n], Ax);
        matrixB.MultMv(&EigVec[i*n], Bx);
        matrixB.MultMv(&EigVec[(i+1)*n], Bx1);
        axpy(n, -EigValR[i], Bx, 1, Ax, 1);
        axpy(n, EigValI[i], Bx1, 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        matrixA.MultMv(&EigVec[(i+1)*n], Ax);
        axpy(n, -EigValI[i], Bx, 1, Ax, 1);
        axpy(n, -EigValR[i], Bx1, 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(EigValR[i], EigValI[i]);
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] Bx1;
    delete[] ResNorm;

  }

} // Solution.


#endif // ANSYMSOL_H

