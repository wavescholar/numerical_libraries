/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE AREig.h.
   Example functions AREig. This function returns eigenvalues, EigVal, and
   eigenvectors, EigVec, of real symmetric, real nonsymmetric and complex
   problems in various modes using ARluSymStdEig, ARluNonSymStdEig and
   ARluNonSymGenEig classes.
   The SuperLU package is employed to solve the linear systems that appear
   when the shift-and-invert spectral transformation is being used.

   There are eighteen different versions of AREig, as shown below. The type
   and the meaning of each AREig parameter is briefly described in section
   II. For a complete description of all parameters, see the Appendix of
   "ARPACK++: a c++ implementation of ARPACK eigenvalue package."

   I) How to call AREig:

      AREig is a c++ overloaded function, which means that there are
      several definitions for this single function. Each definition is
      related to a different problem and a different kind of data. All
      twenty six available AREig functions are listed below. They are
      divided according to the problem type (standard or generalized),
      the type of the involved matrices (real or complex), the
      computational mode being used (regular or shift-and-invert) and
      the desired output data (only eigenvalues or eigenvalues and
      eigenvectors).

      1) Real symmetric standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, uplo, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, uplo, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      2) Real symmetric generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, uplo, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, uplo, nev, which,
                          ncv, tol, maxit, resid, AutoShift)

         b) Shift-and-invert, buckling and Cayley modes:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, uplo, InvertMode, sigma, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA, nnzB,
                          B, irowB, pcolB, uplo, InvertMode, sigma, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

      3) Real nonsymmetric standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol,
                          sigma, nev, which, ncv, tol, maxit, resid, AutoShift)

      4) Real nonsymmetric generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA,
                          pcolA, nnzB, B, irowB, pcolB, nev, which, ncv,
                          tol, maxit, resid, AutoShift)

         b) Real shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, sigma, nev, which, ncv,
                          tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA,
                          pcolA, nnzB, B, irowB, pcolB, sigma, nev, which,
                          ncv, tol, maxit, resid, AutoShift)

         c) Complex shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigValR, EigValI, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, part, sigmaR, SigmaI,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigValR, EigValI, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, part, sigmaR, SigmaI,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      5) Complex standard problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, nev,
                          which, ncv, tol, maxit, resid, AutoShift)

         b) Shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnz, A, irow, pcol, sigma,
                          nev, which, ncv, tol, maxit, resid, AutoShift)

      6) Complex generalized problems:

         a) Regular mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, nev, which, ncv, tol, maxit,
                          resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

         b) Shift-and-invert mode:

            i)  To obtain only eigenvalues:

                int AREig(EigVal, n, nnzA, A, irowA, pcolA, nnzB, B,
                          irowB, pcolB, sigma, nev, which, ncv, tol,
                          maxit, resid, AutoShift)

            ii) To obtain eigenvalues and eigenvectors:

                int AREig(EigVal, EigVec, n, nnzA, A, irowA, pcolA,
                          nnzB, B, irowB, pcolB, sigma, nev, which,
                          ncv, tol, maxit, resid, AutoShift)


   II) AREig parameters:

   Each AREig function contains a variable number of parameters, because
   some parameters are not mandatory. Optional parameters should only be
   defined by the used if the default values supplied by the function are
   not suitable.
   In the description of AREig parameters, the following notation is used
   to describe eigenvalue problems:
   a) standard problems:    A*EigVec = Eigvec*EigVal,
   b) generalized problems: A*EigVec = B*EigVec*EigVal.

      1) Compulsory input parameters:

      int n            Dimension of the problem.

      int nnz          Number of nonzero elements in matrix A.

      int nnzA         Same as nnz.

      int nnzB         Number of nonzero elements in matrix B.

      TYPE A[]         Array of nonzero elements in matrix A.
                       TYPE must be one of "float", "double",
                       "arcomplex<float>" or "arcomplex<double>".

      TYPE B[]         Array of nonzero elements in matrix B.
                       TYPE must be one of "float", "double",
                       "arcomplex<float>" or "arcomplex<double>".

      int irow[]       Array of row indices of the nonzero elements in A.

      int irowA[]      Same as irow.

      int irowB[]      Array of row indices of the nonzero elements in B.

      int pcol[]       Array of pointers to the beginning of columns in
                       A and irow. pcol must have n+1 elements and the
                       last element must be nnz.

      int pcolA[]      Same as pcol.

      int pcolB[]      Array of pointers to the beginning of columns in
                       B and irowB. pcol must have n+1 elements and the
                       last element must be nnzB.

      char uplo        A parameter used only if the problem is symmetric.
                       uplo indicates whether the lower uplo = ’L’) or the
                       upper triangular (uplo = ’U’) part of A (and also
                       B, if the problem is a generalized one) is being
                       supplied by the user.

      char InvertMode  Spectral transformation used when solving symmetric
                       generalized problems. To use the shift and invert
                       mode, the user must set InvertMode to 'S'.  Buckling
                       and Cayley modes are represented by 'B' and 'C',
                       respectively.

      TYPE sigma       The shift (when shift-and invert mode is being used).
                       TYPE must be one of "float" or "double" if the
                       problem is nonsymmetric and sigma is real, or one of
                       "arcomplex<float>" or "arcomplex<double>" if the
                       problem is complex.

      FLOAT sigmaR     Real part of the shift if the problem is real but
                       the shift is complex. FLOAT must be one of "float"
                       or "double".

      FLOAT sigmaI     Imaginary part of the shift if the problem is real
                       but the shift is complex. FLOAT must be one of
                       "float" or "double".

      char part        A parameter that characterizes which part (real or
                       imaginary) of the vector y = OP*x will be used by
                       ARPACK++ when the problem is real but the shift is
                       complex. "part" must be set to one of 'R' (real
                       part) or 'I' (imaginary part).

      int  nev         Number of eigenvalues to be computed.


      2) Optional input parameters:

      char* which      A parameter thar specifies which of the Ritz values
                       are to be computed. "which" must be set to one of:
                       LM: to find eigenvalues with largest magnitude;
                       SM: to find eigenvalues with smallest magnitude;
                       LR: to find eigenvalues with largest real part;
                       SR: to find eigenvalues with smallest real part;
                       LI: to find eigenvalues with largest imaginary part;
                       SI: to find eigenvalues with smallest imaginary part.
                       Default: LM.

      int  ncv         Number of Arnoldi vectors generated at each
                       iteration of ARPACK. Default: 2*nev+1.

      FLOAT tol        Stopping criterion (relative accuracy of Ritz
                       values). FLOAT must be one of "float" or "double".
                       Default: machine precision.

      int  maxit       Maximum number of Arnoldi update iterations allowed.
                       Default: 100*nev.

      TYPE* resid      A pointer to an array that contains the initail
                       vector. Default: a random vector.
                       TYPE must be one of "float", "double",
                       "arcomplex<float>" or "arcomplex<double>".

      bool AutoShift   A parameter that indicates if exact shifts for
                       implicit restarting of the Arnoldi method are to
                       be generated internally by ARPACK++ or shifts are
                       being supplied by the user. Default: true (exact
                       shifts are being used).

      3) Output parameters:

      TYPE  EigVal[]   A vector that contains the "converged" eigenvalues
                       when the problem is symmetric or complex. TYPE
                       must be one float, double, "arcomplex<float>" or
                       "arcomplex<double>".

      FLOAT EigValR[]  A vector that contains the real part of the
                       "converged" eigenvalues (when the problem is real).
                       FLOAT must be one of "float" or "double".

      FLOAT EigValI[]  A vector that contains the imaginary part of the
                       "converged" eigenvalues (when the problem is real).
                       FLOAT must be one of "float" or "double".

      TYPE  EigVec[]   A vector that stores all "converged" eigenvectors
                       consecutively. For real problems, complex eigenvectors
                       are given as two consecutive vectors. The first
                       contains the real part of the eigenvector, while
                       the imaginary part is stored in the second vector.
                       TYPE must be one of "float", "double",
                       "arcomplex<float>" or "arcomplex<double>".


   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef AREIG_H
#define AREIG_H


#include "arcomp.h"
#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnz, arcomplex<FLOAT> A[],
          int irow[], int pcol[], int nev, char* which = "LM", int ncv = 0,
          FLOAT tol = 0.0, int maxit = 0, arcomplex<FLOAT>* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                             maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnz, arcomplex<FLOAT> A[], int irow[], int pcol[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                             maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnz, arcomplex<FLOAT> A[],
          int irow[], int pcol[], arcomplex<FLOAT> sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                             tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnz, arcomplex<FLOAT> A[], int irow[], int pcol[],
          arcomplex<FLOAT> sigma, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluCompStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                             tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnzA,
          arcomplex<FLOAT> A[], int irowA[], int pcolA[], int nnzB,
          arcomplex<FLOAT> B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnzA, arcomplex<FLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<FLOAT> B[], int irowB[], int pcolB[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], int n, int nnzA, arcomplex<FLOAT> A[],
          int irowA[], int pcolA[], int nnzB, arcomplex<FLOAT> B[],
          int irowB[], int pcolB[], arcomplex<FLOAT> sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // complex generalized problem, only eigenvalues, shift-and-invert mode.


template <class FLOAT>
int AREig(arcomplex<FLOAT> EigVal[], arcomplex<FLOAT> EigVec[], int n,
          int nnzA, arcomplex<FLOAT> A[], int irowA[], int pcolA[],
          int nnzB, arcomplex<FLOAT> B[], int irowB[], int pcolB[],
          arcomplex<FLOAT> sigma, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          arcomplex<FLOAT>* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<arcomplex<FLOAT> > matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<arcomplex<FLOAT> > matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluCompGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                             ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // complex generalized problem, values and vectors, shift-and-invert mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                               maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnz,
          FLOAT A[], int irow[], int pcol[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n,
          int nnzA, FLOAT A[], int irowA[], int pcolA[],
          int nnzB, FLOAT B[], int irowB[], int pcolB[],
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(double EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(float EigValR[], FLOAT EigValI[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n,
          int nnzA, FLOAT A[], int irowA[], int pcolA[], int nnzB,
          FLOAT B[], int irowB[], int pcolB[], FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, sigma, which,
                               ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // real shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char part, FLOAT sigmaR, FLOAT sigmaI,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, part,
                               sigmaR, sigmaI, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigValR, EigValI);

} // real nonsymmetric generalized problem, only eigenvalues,
  // complex shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigValR[], FLOAT EigValI[], FLOAT EigVec[], int n, int nnzA,
          FLOAT A[], int irowA[], int pcolA[], int nnzB, FLOAT B[],
          int irowB[], int pcolB[], char part, FLOAT sigmaR, FLOAT sigmaI,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluNonSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA);
  ARluNonSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB);

  // Defining the eigenvalue problem.

  ARluNonSymGenEig<FLOAT> prob(nev, matrixA, matrixB, part,
                               sigmaR, sigmaI, which, ncv,
                               tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigValR, EigValI);

} // real nonsymmetric generalized problem, values and vectors,
  // complex shift-and-invert mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnz, FLOAT A[], int irow[],
          int pcol[], char uplo, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                            maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnz, FLOAT A[],
          int irow[], int pcol[], char uplo, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, which, ncv, tol,
                            maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnz, FLOAT A[], int irow[],
          int pcol[], char uplo, FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                            tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric standard problem, only eigenvalues, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnz, FLOAT A[],
          int irow[], int pcol[], char uplo, FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating a matrix in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  // Defining the eigenvalue problem.

  ARluSymStdEig<FLOAT> prob(nev, matrix, sigma, which, ncv,
                            tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric standard problem, values and vectors, shift-and-invert.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnzA, FLOAT A[], int irowA[],
          int pcolA[], int nnzB, FLOAT B[], int irowB[], int pcolB[],
          char uplo, int nev, char* which = "LM", int ncv = 0,
          FLOAT tol = 0.0, int maxit = 0, FLOAT* resid = 0,
          bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                            ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                            ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors, regular mode.


template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnzA, FLOAT A[], int irowA[],
          int pcolA[], int nnzB, FLOAT B[], int irowB[], int pcolB[],
          char uplo, char InvertMode, FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues,
  // shift-and-invert, buckling and Cayley modes.


template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, char InvertMode, FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors,
  // shift-and-invert, buckling and Cayley modes.


#endif // AREIG_H

