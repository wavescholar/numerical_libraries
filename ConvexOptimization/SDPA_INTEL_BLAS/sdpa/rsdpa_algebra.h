/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/*-------------------------------------------------
  rsdpa_algebra.h
-------------------------------------------------*/

#ifndef __rsdpa_algebra_h__
#define __rsdpa_algebra_h__

#include "rsdpa_struct.h"

class rAl
{
public:
  static double getMinEigenValue(rDenseMatrix& aMat,
				 rVector& eigenVec,
				 rVector& workVec);
  static double getMinEigenValue(rBlockDenseMatrix& aMat,
				 rBlockVector& eigenVec,
				 rBlockVector& workVec);

  static bool getInnerProduct(double& ret,
			      rVector& aVec, rVector& bVec);
  static bool getInnerProduct(double& ret,
			      rBlockVector& aVec,
			      rBlockVector& bVec);
  static bool getInnerProduct(double& ret,
			      rDenseMatrix& aMat,
			      rDenseMatrix& bMat);
  static bool getInnerProduct(double& ret,
			      rSparseMatrix& aMat,
			      rDenseMatrix&  bMat);
  static bool getInnerProduct(double& ret,
			      rBlockDenseMatrix& aMat,
			      rBlockDenseMatrix& bMat);
  static bool getInnerProduct(double& ret,
			      rBlockSparseMatrix& aMat,
			      rBlockDenseMatrix&  bMat);

  static bool getCholesky(rDenseMatrix& retMat, rDenseMatrix& aMat);

  static bool getInvLowTriangularMatrix(rDenseMatrix& retMat,
					rDenseMatrix& aMat);
  static bool getCholeskyAndInv(rBlockDenseMatrix& choleskyMat,
				rBlockDenseMatrix& inverseMat,
				rBlockDenseMatrix& aMat);

  static bool getSymmetrize(rDenseMatrix& aMat);
  static bool getSymmetrize(rBlockDenseMatrix& aMat);

  static bool getTranspose(rDenseMatrix& retMat,
			   rDenseMatrix& aMat);
  static bool getTranspose(rBlockDenseMatrix& retMat,
			   rBlockDenseMatrix& aMat);
  
  static int rdpotf2_(char*uplo, int *n, double *a, int *lda, int *info);
  static int rdpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
  static bool choleskyFactorWithAdjust(rDenseMatrix& aMat);
  
  static bool solveSystems(rVector& xVec,
			   rDenseMatrix& aMat, rVector& bVec);
  // solve aMat * xVec = bVec
  // aMat must be Cholesky Factorized.

  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, rDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rSparseMatrix& aMat, rDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, rSparseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       double* scalar = NULL);
  static bool multiply(rVector& retVec,
		       rVector& aVec, double* scalar = NULL);
  static bool multiply(rBlockVector& retVec,
		       rBlockVector& aVec,
		       double* scalar = NULL);
  static bool multiply(rVector& retVec,
		       rDenseMatrix& aMat, rVector& bVec,
		       double* scalar = NULL);
  static bool multiply(rBlockVector& retVec,
		       rBlockDenseMatrix& aMat,
		       rBlockVector& bVec,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       rBlockDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockSparseMatrix& aMat,
		       rBlockDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       rBlockSparseMatrix& bMat,
		       double* scalar = NULL);
  // ret = aMat**T * bMat
  static bool tran_multiply(rDenseMatrix& retMat,
			    rDenseMatrix& aMat, rDenseMatrix& bMat,
			    double* scalar = NULL);
  static bool tran_multiply(rBlockDenseMatrix& retMat,
			    rBlockDenseMatrix& aMat,
			    rBlockDenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = aMat * bMat**T
  static bool multiply_tran(rDenseMatrix& retMat,
			    rDenseMatrix& aMat, rDenseMatrix& bMat,
			    double* scalar = NULL);
  static bool multiply_tran(rBlockDenseMatrix& retMat,
			    rBlockDenseMatrix& aMat,
			    rBlockDenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = a + (*scalar)*b
  static bool plus(rVector& retVec, rVector& aVec,
		   rVector& bVec, double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rSparseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rSparseMatrix& bMat,
		   double* scalar = NULL);
  
  static bool plus(rBlockVector& retVec,
		   rBlockVector& aVec,
		   rBlockVector& bVec, double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockSparseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockSparseMatrix& bMat,
		   double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rVector& retVec, const char eq,
		  rVector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rBlockVector& retVec, const char eq,
		  rBlockVector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '+' '-' b*(*scalar)
  static bool let(rVector& retVec, const char eq,
		  rVector& aVec, const char op,
		  rVector& bVec, double* scalar = NULL);

  // ret = a '+' '-' '*' 't' 'T' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rSparseMatrix& aMat, const char op,
		  rDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rSparseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' 't' 'T' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockSparseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockSparseMatrix& bMat, double* scalar = NULL);


  // ret = aMat '*' '/' bVec
  static bool let(rVector& rVec, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rVector& bVec);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rVector& aVec, const char op,
		  rVector& bVec);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rDenseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rSparseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rSparseMatrix& aMat, const char op,
		  rDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockVector& aVec, const char op,
		  rBlockVector& bVec);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockSparseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockSparseMatrix& bMat);

};

#endif // __rsdpa_algebra_h__
