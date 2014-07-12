/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SymSol.h
   Template functions that exemplify how to print information 
   about symmetric standard eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SYMSOL_H
#define SYMSOL_H

#include <math.h>
#include "blas1c.h"
#include "matprod.h"
#include "arssym.h"

template<class MATRIX, class FLOAT>
void Solution(MATRIX &A, ARSymStdEig<FLOAT, MATRIX> &Prob)
/*
  Prints eigenvalues and eigenvectors of symmetric eigen-problems
  on standard "cout" stream.
*/

{

  int   i, n, nconv, mode;
  FLOAT *Ax;
  FLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  cout << endl << endl << "Testing ARPACK++ class ARSymStdEig \n";
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << endl;
  switch (mode) {
  case 1:
    cout << "Regular mode" << endl << endl;
    break;
  case 3: 
    cout << "Shift and invert mode" << endl << endl;
  }

  cout << "Dimension of the system            : " << n              << endl;
  cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev()  << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
  cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv()  << endl;
  cout << "Number of iterations taken         : " << Prob.GetIter() << endl;
  cout << endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    cout << "Eigenvalues:" << endl;
    for (i=0; i<nconv; i++) {
      cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << endl;
    }
    cout << endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv];

    for (i=0; i<nconv; i++) {
      A.MultMv(Prob.RawEigenvector(i),Ax);
      axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


#endif // SYMSOL_H
