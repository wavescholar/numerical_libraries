/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE BNMatrxE.h
   Function template for the tridiagonal matrix with 4 on the main
   diagonal and 1 on the subdiagonal and superdiagonal.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef BNMATRXE_H
#define BNMATRXE_H

template<class FLOAT, class INT>
void BandMatrixE(INT n, INT& nL, INT& nU, FLOAT* &A)

{

  INT   j;

  // Defining the upper and lower bandwidth.

  nL  = 1;
  nU  = 1;

  // Creating output vector A.

  A   = new FLOAT[3*n];
  for (j=1; j<(3*n); j+=3) {
    if (j-1)    A[j-1] = 1.0;
                A[j]   = 4.0;
    if (n-j-1)  A[j+1] = 1.0;
  }
    
} //  BandMatrixE.

#endif // BNMATRXE_H

