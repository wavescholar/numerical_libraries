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

#include "rsdpa_algebra.h"
#include "rsdpa_dpotrf.h"

double rAl::getMinEigenValue(rDenseMatrix& aMat,
			     rVector& eigenVec,
			     rVector& workVec)
{
  // aMat is rewritten.
  // aMat must be symmetric.
  // eigenVec is the space of eigen values
  // and needs memory of length aMat.nRow 
  // workVec is temporary space and needs
  // 3*aMat.nRow-1 length memory.
  int N = aMat.nRow;
  int LWORK, info;
  double min_eigen;
  switch (aMat.De_Di) {
  case rDenseMatrix::DENSE:
    LWORK = 3*N-1;
    // "N" means that we need not eigen vectors
    // "L" means that we refer only lower triangular.
    dsyev("NonVectors","Lower",&N,aMat.de_ele,&N,
	   eigenVec.ele,workVec.ele,&LWORK,&info);
    if (info!=0) {
      if (info < 0) {
	rMessage("getMinEigenValue:: info is mistaken " << info);
      } else {
	rMessage("getMinEigenValue:: cannot decomposition");
      }
      exit(0);
      return 0.0;
    }
    return eigenVec.ele[0];
    // Eigen values are sorted by ascending order.
    break;
  case rDenseMatrix::DIAGONAL:
    min_eigen = aMat.di_ele[0];
    eigenVec.ele[0] = aMat.di_ele[0];
    for (int i=1; i<N; ++i) {
      eigenVec.ele[i] = aMat.di_ele[i];
      if (aMat.di_ele[i] < min_eigen) {
	min_eigen = aMat.di_ele[i];
      }
    }
    return min_eigen;
    break;
  }
  return 0.0;
}

double rAl::getMinEigenValue(rBlockDenseMatrix& aMat,
			     rBlockVector& eigenVec,
			     rBlockVector& workVec)
{
  int N = aMat.ele[0].nRow;
  if (eigenVec.ele[0].nDim != N) {
    rError("getMinEigenValue:: different memory size");
  }
  if (workVec.ele[0].nDim != 3*N-1) {
    rError("getMinEigenValue:: different memory size");
  }
  double min_eigen = getMinEigenValue(aMat.ele[0],
				      eigenVec.ele[0],
				      workVec.ele[0]);
  
  for (int l=1; l<aMat.nBlock; ++l) {
    N = aMat.ele[l].nRow;
    if (eigenVec.ele[l].nDim != N) {
      rError("getMinEigenValue:: different memory size");
    }
    // rMessage(" workVec.ele[l].nDim = " << workVec.ele[l].nDim);
    // rMessage(" 3*N-1 = " << 3*N-1);
    if (workVec.ele[l].nDim != 3*N-1) {
      rError("getMinEigenValue:: different memory size");
    }
    double tmp_eigen = getMinEigenValue(aMat.ele[l],
					eigenVec.ele[l],
					workVec.ele[l]);
    if (tmp_eigen < min_eigen) {
      min_eigen = tmp_eigen;
    }
  } // end of for
  return min_eigen;
}

bool rAl::getInnerProduct(double& ret, rVector& aVec, rVector& bVec)
{
   int N = aVec.nDim;
  if (N != bVec.nDim) {
    rError("getInnerProduct:: different memory size");
  }
  ret = ddot(&N,aVec.ele,&IONE,bVec.ele,&IONE);
  ret = ddot(&N,aVec.ele,&IONE,bVec.ele,&IONE);
  return _SUCCESS;
}

bool rAl::getInnerProduct(double& ret,
			  rBlockVector& aVec, rBlockVector& bVec)
{
  if (aVec.nBlock != bVec.nBlock) {
    rError("getInnerProduct:: different memory size");
  }
  bool total_judge = _SUCCESS;
  ret = 0.0;
  double tmp_ret;
  for (int l=0; l<aVec.nBlock; ++l) {
    bool judge = getInnerProduct(tmp_ret,aVec.ele[l],bVec.ele[l]);
    ret += tmp_ret;
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::getInnerProduct(double& ret,
			  rDenseMatrix& aMat, rDenseMatrix& bMat)
{
  if (aMat.nRow!=bMat.nRow || aMat.nCol!=bMat.nCol) {
    rError("getInnerProduct:: different memory size");
  }
  int length;
  switch (aMat.De_Di) {
  case rDenseMatrix::DENSE:
    length = aMat.nRow*aMat.nCol;
    ret = ddot(&length,aMat.de_ele,&IONE,bMat.de_ele,&IONE);
    break;
  case rDenseMatrix::DIAGONAL:
    ret = ddot(&aMat.nCol,aMat.di_ele,&IONE,bMat.di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}

bool rAl::getInnerProduct(double& ret,
			  rSparseMatrix& aMat, rDenseMatrix& bMat)
{
  if (aMat.nRow!=bMat.nRow || aMat.nCol!=bMat.nCol) {
    rError("getInnerProduct:: different memory size");
  }
  int length;
  int amari,shou;
  int index=0, counter=0;
  
  switch(aMat.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    // Attension: in SPARSE case, only half elements
    // are stored. And bMat must be DENSE case.
    ret = 0.0;
    // rMessage("aMat.NonZeroCount == " << aMat.NonZeroCount);
    #if 0
    for (index=0; index<aMat.NonZeroCount; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      double value = aMat.sp_ele      [index];
      // rMessage("i=" << i << "  j=" << j);
      if (i==j) {
	ret+= value*bMat.de_ele[i+bMat.nRow*j];
      } else {
	ret+= value*(bMat.de_ele[i+bMat.nRow*j]
		     + bMat.de_ele[j+bMat.nRow*i]);

      }
    }
    #else
    amari = aMat.NonZeroCount % 4;
    shou = aMat.NonZeroCount / 4;
    for (index=0; index<amari; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      double value = aMat.sp_ele      [index];
      // rMessage("i=" << i << "  j=" << j);
      if (i==j) {
	ret+= value*bMat.de_ele[i+bMat.nRow*j];
      } else {
	ret+= value*(bMat.de_ele[i+bMat.nRow*j]
		     + bMat.de_ele[j+bMat.nRow*i]);

      }
    }

    for (index=amari,counter = 0;
	 counter < shou ; ++counter, index+=4) {
      int        i1 = aMat.row_index   [index];
      int        j1 = aMat.column_index[index];
      double value1 = aMat.sp_ele      [index];
      double ret1;
      // rMessage("i=" << i << "  j=" << j);
      if (i1==j1) {
	ret1= value1*bMat.de_ele[i1+bMat.nRow*j1];
      } else {
	ret1= value1*(bMat.de_ele[i1+bMat.nRow*j1]
		     + bMat.de_ele[j1+bMat.nRow*i1]);

      }
      int        i2 = aMat.row_index   [index+1];
      int        j2 = aMat.column_index[index+1];
      double value2 = aMat.sp_ele      [index+1];
      double ret2;
      // rMessage("i=" << i << "  j=" << j);
      if (i2==j2) {
	ret2= value2*bMat.de_ele[i2+bMat.nRow*j2];
      } else {
	ret2= value2*(bMat.de_ele[i2+bMat.nRow*j2]
		     + bMat.de_ele[j2+bMat.nRow*i2]);

      }
      int        i3 = aMat.row_index   [index+2];
      int        j3 = aMat.column_index[index+2];
      double value3 = aMat.sp_ele      [index+2];
      double ret3;
      // rMessage("i=" << i << "  j=" << j);
      if (i3==j3) {
	ret3= value3*bMat.de_ele[i3+bMat.nRow*j3];
      } else {
	ret3= value3*(bMat.de_ele[i3+bMat.nRow*j3]
		     + bMat.de_ele[j3+bMat.nRow*i3]);

      }
      int        i4 = aMat.row_index   [index+3];
      int        j4 = aMat.column_index[index+3];
      double value4 = aMat.sp_ele      [index+3];
      double ret4;
      // rMessage("i=" << i << "  j=" << j);
      if (i4==j4) {
	ret4= value4*bMat.de_ele[i4+bMat.nRow*j4];
      } else {
	ret4= value4*(bMat.de_ele[i4+bMat.nRow*j4]
		     + bMat.de_ele[j4+bMat.nRow*i4]);

      }
      // ret += ret1;
      // ret += ret2;
      // ret += ret3;
      // ret += ret4;
      ret += (ret1+ret2+ret3+ret4);
    }
    #endif
    break;
  case rSparseMatrix::DENSE:
    length = aMat.nRow*aMat.nCol;
    ret = ddot(&length,aMat.de_ele,&IONE,bMat.de_ele,&IONE);
    break;
  case rSparseMatrix::DIAGONAL:
    ret = ddot(&aMat.nCol,aMat.di_ele,&IONE,bMat.di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}

bool rAl::getInnerProduct(double& ret,
			  rBlockDenseMatrix& aMat,
			  rBlockDenseMatrix& bMat)
{
  if (aMat.nBlock != bMat.nBlock) {
    rError("getInnerProduct:: different memory size");
  }
  bool total_judge = _SUCCESS;
  ret = 0.0;
  double tmp_ret;
  for (int l=0; l<aMat.nBlock; ++l) {
    bool judge = getInnerProduct(tmp_ret,aMat.ele[l],bMat.ele[l]);
    ret += tmp_ret;
    if (judge == FAILURE) {
      rMessage(" something failed");
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::getInnerProduct(double& ret,
			  rBlockSparseMatrix& aMat,
			  rBlockDenseMatrix& bMat)
{
  if (aMat.nBlock != bMat.nBlock) {
    rError("getInnerProduct:: different memory size");
  }
  bool total_judge = _SUCCESS;
  ret = 0.0;
  double tmp_ret;
  for (int l=0; l<aMat.nBlock; ++l) {
    bool judge = getInnerProduct(tmp_ret,aMat.ele[l],bMat.ele[l]);
    ret += tmp_ret;
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::getCholesky(rDenseMatrix& retMat,rDenseMatrix& aMat)
{
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=aMat.nCol
      || retMat.De_Di!=aMat.De_Di) {
    rError("getCholesky:: different memory size");
  }
  int length,info,shou,amari,count;
  int j=0;
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    length = retMat.nRow * retMat.nCol;
    dcopy(&length,aMat.de_ele,&IONE,retMat.de_ele,&IONE);
    #if 1
    dpotrf("Lower",&retMat.nRow,retMat.de_ele,&retMat.nRow,&info);
    #else
    info = choleskyFactorWithAdjust(retMat);
    #endif
    if (info!=0) {
      rMessage("cannot cholesky decomposition");
      rMessage("Could you try with smaller gammaStar?");
      return FAILURE;
    }
    // Make matrix as lower triangular matrix
    #if 0
    for (int j=0; j<retMat.nCol; ++j) {
      for (int i=0; i<j; ++i) {
	retMat.de_ele[i+retMat.nCol*j] = 0.0;
      }
    }
    #else
    for (j=0; j<retMat.nCol; ++j) {
      shou = j /4;
      amari = j%4;
      for (int i=0; i<amari; ++i) {
	retMat.de_ele[i+retMat.nCol*j] = 0.0;
      }
	  int count;
      for (int i=amari,count=0; count < shou; ++count, i+=4) {
	retMat.de_ele[i+retMat.nCol*j] = 0.0;
	retMat.de_ele[i+1+retMat.nCol*j] = 0.0;
	retMat.de_ele[i+2+retMat.nCol*j] = 0.0;
	retMat.de_ele[i+3+retMat.nCol*j] = 0.0;
      }
    }
    #endif
    break;
  case rDenseMatrix::DIAGONAL:
    #if 0
    for (int j=0; j<aMat.nCol; ++j) {
      retMat.di_ele[j] = sqrt(aMat.di_ele[j]);
    }
    #else
    shou = aMat.nCol/4;
    amari = aMat.nCol%4;
    for (j=0; j<amari; ++j) {
      retMat.di_ele[j] = sqrt(aMat.di_ele[j]);
    }
    for (j=amari,count=0; count<shou ; ++count, j+=4) {
      retMat.di_ele[j] = sqrt(aMat.di_ele[j]);
      retMat.di_ele[j+1] = sqrt(aMat.di_ele[j+1]);
      retMat.di_ele[j+2] = sqrt(aMat.di_ele[j+2]);
      retMat.di_ele[j+3] = sqrt(aMat.di_ele[j+3]);
    }
    #endif
    break;
  }
  return _SUCCESS;
}

bool rAl::getInvLowTriangularMatrix(rDenseMatrix& retMat,
				    rDenseMatrix& aMat)
{
  // Make inverse with refference only to lower triangular.
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=aMat.nCol
      || retMat.De_Di!=aMat.De_Di) {
    rError("getCholesky:: different memory size");
  }
  int shou,amari;
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    retMat.setIdentity();
    dtrsm("Left","Lower","NoTraspose","NonUnitDiagonal",
	   &aMat.nRow, &aMat.nCol, &DONE, aMat.de_ele,
	   &aMat.nRow, retMat.de_ele, &retMat.nRow);
    break;
  case rDenseMatrix::DIAGONAL:
    #if 0
    for (int j=0; j<aMat.nCol; ++j) {
      retMat.di_ele[j] = 1.0/aMat.di_ele[j];
    }
    #else
    shou = aMat.nCol/4;
    amari = aMat.nCol%4;
    for (int j=0; j<amari; ++j) {
      retMat.di_ele[j] = 1.0/aMat.di_ele[j];
    }
	int counter;
    for (int j=amari,counter=0; counter<shou; ++counter, j+=4) {
      retMat.di_ele[j] = 1.0/aMat.di_ele[j];
      retMat.di_ele[j+1] = 1.0/aMat.di_ele[j+1];
      retMat.di_ele[j+2] = 1.0/aMat.di_ele[j+2];
      retMat.di_ele[j+3] = 1.0/aMat.di_ele[j+3];
    }
    #endif
    break;
  }
  return _SUCCESS;
}

bool rAl::getCholeskyAndInv(rBlockDenseMatrix& choleskyMat,
			    rBlockDenseMatrix& inverseMat,
			    rBlockDenseMatrix& aMat)
{
  if (choleskyMat.nBlock!=aMat.nBlock
      || inverseMat.nBlock!=aMat.nBlock) {
    rError("getCholeskyAndInv:: different memory size");
  }
  for (int l=0; l<aMat.nBlock; ++l) {
    if (getCholesky(choleskyMat.ele[l],aMat.ele[l]) == false) {
      return false;
    }
    getInvLowTriangularMatrix(inverseMat.ele[l],choleskyMat.ele[l]);
  }
  return _SUCCESS;
}

bool rAl::getSymmetrize(rDenseMatrix& aMat)
{
	int index = 0;
  switch (aMat.De_Di) {
  case rDenseMatrix::DENSE:
    if (aMat.nRow != aMat.nCol) {
      rError("getSymmetrize:: different memory size");
    }
    for (index = 0; index<aMat.nRow-1; ++index) {
      int index1 = index+index*aMat.nRow + 1;
      int index2 = index+(index+1)*aMat.nRow;
      int length = aMat.nRow - 1 - index;
      // aMat.de_ele[index1] += aMat.de_ele[index2]
      daxpy(&length,&DONE,&aMat.de_ele[index2],&aMat.nRow,
	     &aMat.de_ele[index1],&IONE);
      // aMat.de_ele[index1] /= 2.0
      double half = 0.5;
      dscal(&length,&half,&aMat.de_ele[index1],&IONE);
      // aMat.de_ele[index2] = aMat.de_ele[index1]
      dcopy(&length,&aMat.de_ele[index1],&IONE,
	     &aMat.de_ele[index2],&aMat.nRow);
    }
    break;
  case rDenseMatrix::DIAGONAL:
    // Nothing needs.
    break;
  }
  return _SUCCESS;
}

bool rAl::getSymmetrize(rBlockDenseMatrix& aMat)
{
  bool total_judge = _SUCCESS;
  for (int l=0; l<aMat.nBlock; ++l) {
    bool judge = getSymmetrize(aMat.ele[l]);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::getTranspose(rDenseMatrix& retMat,
		       rDenseMatrix& aMat)
{
  if (aMat.nRow != aMat.nCol) {
    rError("getTranspose:: different memory size");
    // Of course, a non-symmetric matrix has
    // its transposed matrix,
    // but in this algorithm we have to make
    // transposed matrix only when symmetric matrix.
  }
  retMat.copyFrom(aMat);
  int i=0;
  switch (aMat.De_Di) {
  case rDenseMatrix::DENSE:
    #if 0
    for (int i=0; i<aMat.nRow; ++i) {
      for (int j=0; j<=i; ++j) {
	int index1 = i+aMat.nCol*j;
	int index2 = j+aMat.nCol*i;
	retMat.de_ele[index1] = aMat.de_ele[index2];
	retMat.de_ele[index2] = aMat.de_ele[index1];
      }
    }
    #else
    for (i=0; i<aMat.nRow; ++i) {
      int shou = (i+1)/4;
      int amari = (i+1)/4;
      for (int j=0; j<amari; ++j) {
	int index1 = i+aMat.nCol*j;
	int index2 = j+aMat.nCol*i;
	retMat.de_ele[index1] = aMat.de_ele[index2];
	retMat.de_ele[index2] = aMat.de_ele[index1];
      }
	  int counter;
      for (int j=amari,counter =0 ; counter < shou;
	   ++counter, j+=4) {
	int index1 = i+aMat.nCol*j;
	int index_1 = j+aMat.nCol*i;
	retMat.de_ele[index1] = aMat.de_ele[index_1];
	retMat.de_ele[index_1] = aMat.de_ele[index1];
	int index2 = i+aMat.nCol*(j+1);
	int index_2 = (j+1)+aMat.nCol*i;
	retMat.de_ele[index2] = aMat.de_ele[index_2];
	retMat.de_ele[index_2] = aMat.de_ele[index2];
	int index3 = i+aMat.nCol*(j+2);
	int index_3 = (j+2)+aMat.nCol*i;
	retMat.de_ele[index3] = aMat.de_ele[index_3];
	retMat.de_ele[index_3] = aMat.de_ele[index3];
	int index4 = i+aMat.nCol*(j+3);
	int index_4 = (j+3)+aMat.nCol*i;
	retMat.de_ele[index4] = aMat.de_ele[index_4];
	retMat.de_ele[index_4] = aMat.de_ele[index4];
      }
    }
    #endif
    break;
  case rDenseMatrix::DIAGONAL:
    // Nothing needs.
    break;
  }
  return _SUCCESS;
}

bool rAl::getTranspose(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat)
{
  if (retMat.nBlock!=aMat.nBlock) {
    rError("getTranspose:: different memory size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<aMat.nBlock; ++l) {
    bool judge = getTranspose(retMat.ele[l],aMat.ele[l]);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

int rAl::rdpotf2_(char*uplo, int *n, double *a, int *lda, int *info)
{
  int nRow = *lda;
  for (int j = 0; j <*n; ++j) {
    double ajj = a[j +nRow*j]
      - ddot(&j, &a[j], lda, &a[j], lda);

    // Here is point.(start)
    if (ajj <= (float)-1.0e-6) {
      a[j + j * nRow] = ajj;
      *info = j+1;
      return 0;
    }
    if (ajj <= (float)1.0e-14) { 
      ajj = 1e100;
      a[j + j * nRow] = ajj;
    } else {
      ajj = sqrt(ajj);
      a[j + j * nRow] = ajj;
    }
    // Here is point.(end)

    if (j < *n-1) {
      int i = *n-1 - j;
      dgemv("No transpose", &i, &j, &DMONE, &a[j + 1],
	     lda, &a[j], lda, &DONE,
	     &a[(j + 1)+nRow*j], &IONE);
      double d1 = 1.0 / ajj;
      dscal(&i, &d1, &a[(j + 1)+nRow*j], &IONE);
    }
  }
  return 0;
}
extern integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len);
int rAl::rdpotrf_(char *uplo, int *n, double *a, int *lda, int *info)
{
  // This funciton makes Cholesky factorization
  // in only case Lower Triangular.
  // That is, A will be L*L**T, not U**T*U.
  int nRow = *lda;
  *info = 0;

  //bbcrevisit issue with num parameters
//MKL_INT     ilaenv_( const MKL_INT *ispec, 
//const char *name, 
//const char *opts, 
//const MKL_INT *n1, 
//const MKL_INT *n2, 
//const MKL_INT *n3, 
//const MKL_INT *n4 );
//  INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
//    ISPEC   (input) INTEGER
//*          Specifies the parameter to be returned as the value of
//*          ILAENV.
//*          = 1: the optimal blocksize; if this value is 1, an unblocked
//*               algorithm will give the best performance.
//*          = 2: the minimum block size for which the block routine
//*               should be used; if the usable block size is less than
//*               this value, an unblocked routine should be used.
//*          = 3: the crossover point (in a block routine, for N less
//*               than this value, an unblocked routine should be used)
//*          = 4: the number of shifts, used in the nonsymmetric
//*               eigenvalue routines (DEPRECATED)
//*          = 5: the minimum column dimension for blocking to be used;
//*               rectangular blocks must have dimension at least k by m,
//*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
//*          = 6: the crossover point for the SVD (when reducing an m by n
//*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
//*               this value, a QR factorization is used first to reduce
//*               the matrix to a triangular form.)
//*          = 7: the number of processors
//*          = 8: the crossover point for the multishift QR method
//*               for nonsymmetric eigenvalue problems (DEPRECATED)
//*          = 9: maximum size of the subproblems at the bottom of the
//*               computation tree in the divide-and-conquer algorithm
//*               (used by xGELSD and xGESDD)
//*          =10: ieee NaN arithmetic can be trusted not to trap
//*          =11: infinity arithmetic can be trusted not to trap
//*          12 <= ISPEC <= 16:
//*               xHSEQR or one of its subroutines,
//*               see IPARMQ for detailed explanation
//*
//*  NAME    (input) CHARACTER*(*)
//*          The name of the calling subroutine, in either upper case or
//*          lower case.
//*
//*  OPTS    (input) CHARACTER*(*)
//*          The character options to the subroutine NAME, concatenated
//*          into a single character string.  For example, UPLO = 'U',
//*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
//*          be specified as OPTS = 'UTN'.
//*
//*  N1      (input) INTEGER
//*  N2      (input) INTEGER
//*  N3      (input) INTEGER
//*  N4      (input) INTEGER
//*          Problem dimensions for the subroutine NAME; these may not all
//*          be required.
//*
//* (ILAENV) (output) INTEGER
//*          >= 0: the value of the parameter specified by ISPEC
//*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
//*


//ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
  int nb = ilaenv(&IONE, "DPOTRF", "L", n,  &IMONE,&IONE, &IMONE);//bbcrevisit this horror , 6, 1);
  if (nb <= 1 || nb >= *n) {
    // Here is point.
    rdpotf2_(uplo, n, a, lda, info);
  } else {

    for (int j = 0; j < *n; j += nb) {
      int jb = min(nb,*n- j);
      dsyrk("Lower", "No transpose", &jb,
	     &j, &DMONE, &a[j], lda,
	     &DONE, &a[j+nRow*j], lda);
      // Here is point.
      rdpotf2_("Lower", &jb, &a[j+nRow*j], lda, info);
      if (*info != 0) {
	*info = *info + j - 1;
	return 0;
      }
      if (j + jb <= *n-1) {
	int i = *n - j - jb;
	dgemm("No transpose", "Transpose", &i, &jb,
	       &j, &DMONE, &a[j + jb], lda, &a[j], lda, 
	       &DONE, &a[(j + jb)+nRow*j], lda);
	dtrsm("Right", "Lower", "Transpose", "Non-unit",
	       &i, &jb, &DONE, &a[j+nRow*j], lda,
	       &a[(j + jb)+nRow*j], lda);
      }
    }
  }
  return 0;
}


bool rAl::choleskyFactorWithAdjust(rDenseMatrix& aMat)
{
  int info=0;
#if 1
  // aMat.display();
  rTimeStart(START1);
  info = rATL_dpotrfL(aMat.nRow, aMat.de_ele,aMat.nRow);
  rTimeEnd(END1);
  // rMessage("Schur colesky  ::"  << rTimeCal(START1,END1));
  // aMat.display();
#elif 1
  dpotrf_("Lower",&aMat.nRow,aMat.de_ele,&aMat.nRow,&info);
#else
  rdpotrf_("Lower",&aMat.nRow,aMat.de_ele,&aMat.nRow,&info);
#endif
  if (info < 0) {
    rMessage("cholesky argument is wrong " << -info);
  } else if (info > 0) {
    rMessage("cholesky miss condition :: not positive definite"
	     << " :: info = " << info);

    return FAILURE;
  }
  return _SUCCESS;
#if 0
  double ZERO_DETECT = 1.0e-3;
  double NONZERO = 1.0e-7;
  // no idea version
  // if Cholesky factorization failed, then exit soon.
  int info = 1; // info == 0 means success
  int start = 0;
  while (start<aMat.nRow) {
    int N = aMat.nRow - start;
    dpotf2_("Lower",&N,&aMat.de_ele[start+start*aMat.nRow],
	    &aMat.nRow,&info);
    if (info <=0) {
      // rMessage("Cholesky is very nice");
      break;
    }
    start += (info-1); // next target
    double wrong = aMat.de_ele[start+start*aMat.nRow];
    if (wrong < -ZERO_DETECT) {
      rMessage("cholesky adjust position " << start);
      rMessage("cannot cholesky decomposition"
	       " with adjust " << wrong);
      return FAILURE;
    }
    aMat.de_ele[start+start*aMat.nRow] = NONZERO;
    if (start<aMat.nRow-1) {
      // improve the right down element of 0
      for (int j=1; j<=aMat.nRow-1-start; ++j) {
	double& migi  = aMat.de_ele[start+(start+j)*aMat.nRow];
	double& shita = aMat.de_ele[(start+j)+start*aMat.nRow];
	double& mishi = aMat.de_ele[(start+j)+(start+j)*aMat.nRow];
	// rMessage(" mishi = " << mishi);
	if (mishi < NONZERO) {
	  // rMessage(" mishi < NONZERO ");
	  mishi = NONZERO;
	  migi  = NONZERO * 0.1;
	  shita = NONZERO * 0.1;
	} else if (migi*shita > NONZERO*mishi) {
	  // rMessage(" migi*migi > NONZERO*mishi ");
	  migi  = sqrt(NONZERO*mishi) * 0.99;
	  shita = sqrt(NONZERO*mishi) * 0.99;
	}
      }
    }
    rMessage("cholesky adjust position " << start);
  }
  if (info < 0) {
    rError("argument is something wrong " << info);
  }
  return _SUCCESS;
#endif
}

bool rAl::solveSystems(rVector& xVec,
		       rDenseMatrix& aMat, rVector& bVec)
{
  // aMat must have done Cholesky factorized.
  if (aMat.nCol!=xVec.nDim || aMat.nRow!=bVec.nDim
      || aMat.nRow!=aMat.nCol) {
    rError("solveSystems:: different memory size");
  }
  if (aMat.De_Di!=rDenseMatrix::DENSE) {
    rError("solveSystems:: matrix type must be DENSE");
  }
  xVec.copyFrom(bVec);
  dtrsv("Lower", "NoTranspose", "NonUnit",
	 &aMat.nRow, aMat.de_ele, &aMat.nCol, xVec.ele,&IONE);
  dtrsv("Lower", "Transpose", "NonUnit",
	 &aMat.nRow, aMat.de_ele, &aMat.nCol, xVec.ele,&IONE);
  return _SUCCESS;
}

bool rAl::multiply(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar)
{
  if (retMat.nRow!=aMat.nRow || aMat.nCol!=bMat.nRow
      || bMat.nCol!=retMat.nCol
      || retMat.De_Di!=aMat.De_Di || retMat.De_Di!=bMat.De_Di) {
    rError("multiply :: different matrix size");
  }
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    if (scalar==NULL) {
      scalar = &DONE;
      // attension::scalar is loval variable.
    }
    dgemm("NoTranspose","NoTranspose",
	   &retMat.nRow,&retMat.nCol,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nRow,bMat.de_ele,&bMat.nRow,
	   &DZERO,retMat.de_ele,&retMat.nRow);
    break;
  case rDenseMatrix::DIAGONAL:
    if (scalar==NULL) {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int  j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    }
    break;
  }
  return _SUCCESS;
}

bool rAl::multiply(rDenseMatrix& retMat,
		   rSparseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar)
{
  if (retMat.nRow!=aMat.nRow || aMat.nCol!=bMat.nRow
      || bMat.nCol!=retMat.nCol) {
    rError("multiply :: different matrix size");
  }
  retMat.setZero();
  switch (aMat.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| bMat.De_Di!=rDenseMatrix::DENSE) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      for (int index=0; index<aMat.NonZeroCount; ++index) {
	int        i = aMat.row_index   [index];
	int        j = aMat.column_index[index];
	double value = aMat.sp_ele      [index];
	if (i!=j) {
	  #define MULTIPLY_NON_ATLAS 1
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[i+retMat.nRow*t] +=
	      value * bMat.de_ele[j+bMat.nRow*t];
	    retMat.de_ele[j+retMat.nRow*t] +=
	      value * bMat.de_ele[i+bMat.nRow*t];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[j],&bMat.nRow,
		 &retMat.de_ele[i],&retMat.nRow);
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[i],&bMat.nRow,
		 &retMat.de_ele[j],&retMat.nRow);
	  #endif
	} else {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[j+retMat.nRow*t] +=
	      value * bMat.de_ele[j+bMat.nRow*t];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[j],&bMat.nRow,
		 &retMat.de_ele[j],&retMat.nRow);
	  #endif
	}
      } // end of 'for index'
    } else { // scalar!=NULL
      for (int index=0; index<aMat.NonZeroCount; ++index) {
	int        i = aMat.row_index   [index];
	int        j = aMat.column_index[index];
	double value = aMat.sp_ele      [index] * (*scalar);
	if (i!=j) {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[i+retMat.nRow*t] +=
	      value * bMat.de_ele[j+bMat.nRow*t];
	    retMat.de_ele[j+retMat.nRow*t] +=
	      value * bMat.de_ele[i+bMat.nRow*t];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[j],&bMat.nRow,
		 &retMat.de_ele[i],&retMat.nRow);
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[i],&bMat.nRow,
		 &retMat.de_ele[j],&retMat.nRow);
	  #endif
	} else {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[j+retMat.nRow*t] +=
	      value * bMat.de_ele[j+bMat.nRow*t];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&bMat.de_ele[j],&bMat.nRow,
		 &retMat.de_ele[j],&retMat.nRow);
	  #endif
	}
      } // end of 'for index'
    } // end of 'if (scalar==NULL)
    break;
  case rSparseMatrix::DENSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| bMat.De_Di!=rDenseMatrix::DENSE) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      scalar = &DONE;
      // attension:: scalar is local variable.
    }
    dgemm("NoTranspose","NoTranspose",
	   &retMat.nRow,&retMat.nCol,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nRow,bMat.de_ele,&bMat.nRow,
	   &DZERO,retMat.de_ele,&retMat.nRow);
    break;
  case rSparseMatrix::DIAGONAL:
    if (retMat.De_Di!=rDenseMatrix::DIAGONAL
	|| bMat.De_Di!=rDenseMatrix::DIAGONAL) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    }
    break;
  } // end of switch

  return _SUCCESS;
}

bool rAl::multiply(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rSparseMatrix& bMat,
		   double* scalar)
{
  if (retMat.nRow!=aMat.nRow || aMat.nCol!=bMat.nRow
      || bMat.nCol!=retMat.nCol) {
    rError("multiply :: different matrix size");
  }
  retMat.setZero();
  switch (bMat.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    // rMessage("Here will be faster by atlas");
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| aMat.De_Di!=rDenseMatrix::DENSE) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      for (int index=0; index<bMat.NonZeroCount; ++index) {
	int        i = bMat.row_index   [index];
	int        j = bMat.column_index[index];
	double value = bMat.sp_ele      [index];
	if (i!=j) {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[t+retMat.nRow*j] +=
	      value * aMat.de_ele[t+aMat.nRow*i];
	    retMat.de_ele[t+retMat.nRow*i] +=
	      value * aMat.de_ele[t+aMat.nRow*j];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*j],&IONE,
		 &retMat.de_ele[retMat.nRow*i],&IONE);
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*i],&IONE,
		 &retMat.de_ele[retMat.nRow*j],&IONE);
	  #endif
	} else {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[t+retMat.nRow*j] +=
	      value * aMat.de_ele[t+aMat.nRow*j];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*j],&IONE,
		 &retMat.de_ele[retMat.nRow*j],&IONE);
	  #endif
	}
      } // end of 'for index'
    } else { // scalar!=NULL
      for (int index=0; index<bMat.NonZeroCount; ++index) {
	int        i = bMat.row_index   [index];
	int        j = bMat.column_index[index];
	double value = bMat.sp_ele      [index] * (*scalar);
	if (i!=j) {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[t+retMat.nCol*j] +=
	      value * aMat.de_ele[t+bMat.nCol*i];
	    retMat.de_ele[t+retMat.nCol*i] +=
	      value * aMat.de_ele[t+bMat.nCol*j];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*j],&IONE,
		 &retMat.de_ele[retMat.nRow*i],&IONE);
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*i],&IONE,
		 &retMat.de_ele[retMat.nRow*j],&IONE);
	  #endif
	} else {
	  #if MULTIPLY_NON_ATLAS
	  for (int t=0; t<bMat.nCol; ++t) {
	    retMat.de_ele[t+retMat.nCol*j] +=
	      value * aMat.de_ele[t+aMat.nCol*j];
	  }
	  #else
	  daxpy_(&bMat.nCol,&value,&aMat.de_ele[aMat.nRow*j],&IONE,
		 &retMat.de_ele[retMat.nRow*j],&IONE);
	  #endif
	}
      } // end of 'for index'
    } // end of 'if (scalar==NULL)
    break;
  case rSparseMatrix::DENSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| aMat.De_Di!=rDenseMatrix::DENSE) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      scalar = &DONE;
      // attension: scalar is local variable.
    }
    dgemm("NoTranspose","NoTranspose",
	   &retMat.nRow,&retMat.nCol,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nRow,bMat.de_ele,&bMat.nRow,
	   &DZERO,retMat.de_ele,&retMat.nRow);
    break;
  case rSparseMatrix::DIAGONAL:
    if (retMat.De_Di!=rDenseMatrix::DIAGONAL
	|| aMat.De_Di!=rDenseMatrix::DIAGONAL) {
      rError("multiply :: different matrix type");
    }
    if (scalar==NULL) {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int  j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    }
    break;
  } // end of switch

  return _SUCCESS;
}

bool rAl::multiply(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, double* scalar)
{
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=retMat.nCol
      || retMat.De_Di!=aMat.De_Di) {
    rError("multiply :: different matrix size");
  }
  if (scalar==NULL) {
    scalar = &DONE;
  }
  int length;
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    length = retMat.nRow*retMat.nCol;
    dcopy(&length,aMat.de_ele,&IONE,retMat.de_ele,&IONE);
    dscal(&length,scalar,retMat.de_ele,&IONE);
    break;
  case rDenseMatrix::DIAGONAL:
    dcopy(&retMat.nCol,aMat.di_ele,&IONE,retMat.di_ele,&IONE);
    dscal(&retMat.nCol,scalar,retMat.di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}

bool rAl::multiply(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock) {
    rError("multiply:: different memory size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<aMat.nBlock; ++l) {
    bool judge = multiply(retMat.ele[l],aMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply(rVector& retVec,
		   rVector& aVec, double* scalar)
{
  if (retVec.nDim!=aVec.nDim) {
    rError("multiply :: different vector size");
  }
  if (scalar==NULL) {
    scalar = &DONE;
  }
  dcopy(&retVec.nDim,aVec.ele,&IONE,retVec.ele,&IONE);
  dscal(&retVec.nDim,scalar,retVec.ele,&IONE);
  return _SUCCESS;
}

bool rAl::multiply(rBlockVector& retVec,
		   rBlockVector& aVec,
		   double* scalar)
{
  if (retVec.nBlock!=aVec.nBlock) {
    rError("multiply:: different memory size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<aVec.nBlock; ++l) {
    bool judge = multiply(retVec.ele[l],aVec.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply(rVector& retVec,
		   rDenseMatrix& aMat, rVector& bVec,
		   double* scalar)
{
  if (retVec.nDim!=aMat.nRow || aMat.nCol!=bVec.nDim
      || bVec.nDim!=retVec.nDim) {
    rError("multiply :: different matrix size");
  }
  switch (aMat.De_Di) {
  case rDenseMatrix::DENSE:
    if (scalar==NULL) {
      scalar = &DONE;
    }
    dgemv("NoTranspose",&aMat.nRow,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nRow,bVec.ele,&IONE,
	   &DZERO,retVec.ele,&IONE);
    break;
  case rDenseMatrix::DIAGONAL:
    if (scalar==NULL) {
      #if 0
      for (int j=0; j<aMat.nCol; ++j) {
	retVec.ele[j] = aMat.di_ele[j]*bVec.ele[j];
      }
      #else
      int shou = retVec.nDim / 4;
      int amari = retVec.nDim % 4;
      for (int j=0; j<amari; ++j) {
	retVec.ele[j] = aMat.di_ele[j]*bVec.ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retVec.ele[j] = aMat.di_ele[j]*bVec.ele[j];
	retVec.ele[j+1] = aMat.di_ele[j+1]*bVec.ele[j+1];
	retVec.ele[j+2] = aMat.di_ele[j+2]*bVec.ele[j+2];
	retVec.ele[j+3] = aMat.di_ele[j+3]*bVec.ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0; j<aMat.nCol; ++j) {
	retVec.ele[j] = (*scalar) * aMat.di_ele[j]*bVec.ele[j];
      }
      #else
      int shou = retVec.nDim / 4;
      int amari = retVec.nDim % 4;
      for (int j=0; j<amari; ++j) {
	retVec.ele[j] = (*scalar) * aMat.di_ele[j]*bVec.ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retVec.ele[j] = (*scalar) * aMat.di_ele[j]*bVec.ele[j];
	retVec.ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bVec.ele[j+1];
	retVec.ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bVec.ele[j+2];
	retVec.ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bVec.ele[j+3];
      }
      #endif
    }
    break;
  }
  return _SUCCESS;
}

bool rAl::multiply(rBlockVector& retVec,
		   rBlockDenseMatrix& aMat,
		   rBlockVector& bVec,
		   double* scalar)
{
  if (retVec.nBlock!=aMat.nBlock || retVec.nBlock!=bVec.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retVec.nBlock; ++l) {
    bool judge = multiply(retVec.ele[l],aMat.ele[l],
			  bVec.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = multiply(retMat.ele[l],aMat.ele[l],
			  bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply(rBlockDenseMatrix& retMat,
		   rBlockSparseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = multiply(retMat.ele[l],aMat.ele[l],
			  bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockSparseMatrix& bMat,
		   double* scalar )
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = multiply(retMat.ele[l],aMat.ele[l],
			  bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}
  
bool rAl::tran_multiply(rDenseMatrix& retMat,
			rDenseMatrix& aMat, rDenseMatrix& bMat,
			double* scalar)
{
  if (retMat.nRow!=aMat.nCol || aMat.nRow!=bMat.nRow
      || bMat.nCol!=retMat.nCol
      || retMat.De_Di!=aMat.De_Di || retMat.De_Di!=bMat.De_Di) {
    rError("multiply :: different matrix size");
  }
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    if (scalar==NULL) {
      scalar = &DONE;
      // scalar is local variable
    }
    // The Point is the first argument is "Transpose".
    dgemm("Transpose","NoTranspose",
	   &retMat.nRow,&retMat.nCol,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nCol,bMat.de_ele,&bMat.nRow,
	   &DZERO,retMat.de_ele,&retMat.nRow);
    break;
  case rDenseMatrix::DIAGONAL:
    if (scalar==NULL) {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    }
    break;
  }

  return _SUCCESS;
}

bool rAl::tran_multiply(rBlockDenseMatrix& retMat,
			rBlockDenseMatrix& aMat,
			rBlockDenseMatrix& bMat,
			double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = tran_multiply(retMat.ele[l],aMat.ele[l],
			       bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::multiply_tran(rDenseMatrix& retMat,
			rDenseMatrix& aMat, rDenseMatrix& bMat,
			double* scalar)
{
  if (retMat.nRow!=aMat.nRow || aMat.nCol!=bMat.nCol
      || bMat.nRow!=retMat.nRow
      || retMat.De_Di!=aMat.De_Di || retMat.De_Di!=bMat.De_Di) {
    rError("multiply :: different matrix size");
  }
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    if (scalar==NULL) {
      scalar = &DONE;
    }
    // The Point is the first argument is "NoTranspose".
    dgemm("NoTranspose","Transpose",
	   &retMat.nRow,&retMat.nCol,&aMat.nCol,
	   scalar,aMat.de_ele,&aMat.nRow,bMat.de_ele,&bMat.nCol,
	   &DZERO,retMat.de_ele,&retMat.nRow);
    break;
  case rDenseMatrix::DIAGONAL:
    if (scalar==NULL) {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    } else {
      #if 0
      for (int j=0;j<retMat.nCol; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
      #else
      int shou = retMat.nCol / 4;
      int amari = retMat.nCol % 4;
      for (int j=0; j<amari; ++j) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	retMat.di_ele[j] = (*scalar) * aMat.di_ele[j]*bMat.di_ele[j];
	retMat.di_ele[j+1] = (*scalar) * aMat.di_ele[j+1]*bMat.di_ele[j+1];
	retMat.di_ele[j+2] = (*scalar) * aMat.di_ele[j+2]*bMat.di_ele[j+2];
	retMat.di_ele[j+3] = (*scalar) * aMat.di_ele[j+3]*bMat.di_ele[j+3];
      }
      #endif
    }
    break;
  }
  return _SUCCESS;
}

bool rAl::multiply_tran(rBlockDenseMatrix& retMat,
			rBlockDenseMatrix& aMat,
			rBlockDenseMatrix& bMat,
			double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("multiply:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = multiply_tran(retMat.ele[l],aMat.ele[l],
			       bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::plus(rVector& retVec, rVector& aVec,
	       rVector& bVec, double* scalar)
{
  if (retVec.nDim!=aVec.nDim || aVec.nDim!=bVec.nDim) {
    rError("plus :: different matrix size");
  }
  if (scalar==NULL) {
    scalar = &DONE;
  }
  dcopy(&retVec.nDim,aVec.ele,&IONE,retVec.ele,&IONE);
  daxpy(&retVec.nDim,scalar,bVec.ele,&IONE,retVec.ele,&IONE);
  return _SUCCESS;
}

bool rAl::plus(rDenseMatrix& retMat,
	       rDenseMatrix& aMat, rDenseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=aMat.nCol
      || retMat.nRow!=bMat.nRow || retMat.nCol!=bMat.nCol
      || retMat.De_Di!=aMat.De_Di || retMat.De_Di!=bMat.De_Di) {
    rError("plus :: different matrix size");
  }
  if (scalar==NULL) {
      scalar = &DONE;
  }
  int length;
  switch (retMat.De_Di) {
  case rDenseMatrix::DENSE:
    length = retMat.nRow*retMat.nCol;
    dcopy(&length,aMat.de_ele,&IONE,retMat.de_ele,&IONE);
    daxpy(&length,scalar,bMat.de_ele,&IONE,retMat.de_ele,&IONE);
    break;
  case rDenseMatrix::DIAGONAL:
    dcopy(&retMat.nCol,aMat.di_ele,&IONE,retMat.di_ele,&IONE);
    daxpy(&retMat.nCol,scalar,bMat.di_ele,&IONE,
	   retMat.di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}

bool rAl::plus(rDenseMatrix& retMat,
	       rSparseMatrix& aMat, rDenseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=aMat.nCol
      || retMat.nRow!=bMat.nRow || retMat.nCol!=bMat.nCol) {
    rError("plus :: different matrix size");
  }
  // ret = (*scalar) * b
  if (multiply(retMat,bMat,scalar) == FAILURE) {
    return FAILURE;
  }
  int length;
  // ret += a
  int shou,amari;
	int counter;
	int index;
  switch (aMat.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| bMat.De_Di!=rDenseMatrix::DENSE) {
      rError("plus :: different matrix type");
    }
    #if 0
    for (int index=0; index<aMat.NonZeroCount; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      double value = aMat.sp_ele      [index];
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
    #else
    shou = aMat.NonZeroCount / 4;
    amari = aMat.NonZeroCount % 4;
    for (index=0; index<amari; ++index) {
      int        i = aMat.row_index   [index];
      int        j = aMat.column_index[index];
      double value = aMat.sp_ele      [index];
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
    for (index=amari,counter=0;
	 counter<shou; ++counter,index+=4) {
      int        i1 = aMat.row_index   [index];
      int        j1 = aMat.column_index[index];
      double value1 = aMat.sp_ele      [index];
      if (i1!=j1) {
	retMat.de_ele[i1+retMat.nCol*j1] += value1;
	retMat.de_ele[j1+retMat.nCol*i1] += value1;
      } else {
	retMat.de_ele[i1+retMat.nCol*i1] += value1;
      }
      int        i2 = aMat.row_index   [index+1];
      int        j2 = aMat.column_index[index+1];
      double value2 = aMat.sp_ele      [index+1];
      if (i2!=j2) {
	retMat.de_ele[i2+retMat.nCol*j2] += value2;
	retMat.de_ele[j2+retMat.nCol*i2] += value2;
      } else {
	retMat.de_ele[i2+retMat.nCol*i2] += value2;
      }
      int        i3 = aMat.row_index   [index+2];
      int        j3 = aMat.column_index[index+2];
      double value3 = aMat.sp_ele      [index+2];
      if (i3!=j3) {
	retMat.de_ele[i3+retMat.nCol*j3] += value3;
	retMat.de_ele[j3+retMat.nCol*i3] += value3;
      } else {
	retMat.de_ele[i3+retMat.nCol*i3] += value3;
      }
      int        i4 = aMat.row_index   [index+3];
      int        j4 = aMat.column_index[index+3];
      double value4 = aMat.sp_ele      [index+3];
      if (i4!=j4) {
	retMat.de_ele[i4+retMat.nCol*j4] += value4;
	retMat.de_ele[j4+retMat.nCol*i4] += value4;
      } else {
	retMat.de_ele[i4+retMat.nCol*i4] += value4;
      }
    } // end of 'for index'
    #endif
    break;
  case rSparseMatrix::DENSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| bMat.De_Di!=rDenseMatrix::DENSE) {
      rError("plus :: different matrix type");
    }
    length = retMat.nRow*retMat.nCol;
    daxpy(&length,&DONE,aMat.de_ele,&IONE,retMat.de_ele,&IONE);
    break;
  case rSparseMatrix::DIAGONAL:
    if (retMat.De_Di!=rDenseMatrix::DIAGONAL
	|| bMat.De_Di!=rDenseMatrix::DIAGONAL) {
      rError("plus :: different matrix type");
    }
    daxpy(&retMat.nCol,&DONE,aMat.di_ele,&IONE,
	   retMat.di_ele,&IONE);
    break;
  } // end of switch
  return _SUCCESS;
}

bool rAl::plus(rDenseMatrix& retMat,
	       rDenseMatrix& aMat, rSparseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nRow!=aMat.nRow || retMat.nCol!=aMat.nCol
      || retMat.nRow!=bMat.nRow || retMat.nCol!=bMat.nCol) {
    rError("plus :: different matrix size");
  }
  // ret = a
  if (retMat.copyFrom(aMat) == FAILURE) {
    return FAILURE;
  }
  if (scalar==NULL) {
    scalar = &DONE;
  }
  int length,shou,amari,index=0;
  // ret += (*scalar) * b
  switch (bMat.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| aMat.De_Di!=rDenseMatrix::DENSE) {
      rError("plus :: different matrix type");
    }
    #if 0
    for (int index=0; index<bMat.NonZeroCount; ++index) {
      int        i = bMat.row_index   [index];
      int        j = bMat.column_index[index];
      double value = bMat.sp_ele      [index] * (*scalar);
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
    #else
    shou = bMat.NonZeroCount / 4;
    amari = bMat.NonZeroCount % 4;
    for (index=0; index<amari; ++index) {
      int        i = bMat.row_index   [index];
      int        j = bMat.column_index[index];
      double value = bMat.sp_ele      [index] * (*scalar);
      if (i!=j) {
	retMat.de_ele[i+retMat.nCol*j] += value;
	retMat.de_ele[j+retMat.nCol*i] += value;
      } else {
	retMat.de_ele[i+retMat.nCol*i] += value;
      }
    } // end of 'for index'
	int counter;
    for (index=amari,counter=0;
	 counter<shou; ++counter,index+=4) {
      int        i1 = bMat.row_index   [index];
      int        j1 = bMat.column_index[index];
      double value1 = bMat.sp_ele      [index] * (*scalar);
      if (i1!=j1) {
	retMat.de_ele[i1+retMat.nCol*j1] += value1;
	retMat.de_ele[j1+retMat.nCol*i1] += value1;
      } else {
	retMat.de_ele[i1+retMat.nCol*i1] += value1;
      }
      int        i2 = bMat.row_index   [index+1];
      int        j2 = bMat.column_index[index+1];
      double value2 = bMat.sp_ele      [index+1] * (*scalar);
      if (i2!=j2) {
	retMat.de_ele[i2+retMat.nCol*j2] += value2;
	retMat.de_ele[j2+retMat.nCol*i2] += value2;
      } else {
	retMat.de_ele[i2+retMat.nCol*i2] += value2;
      }
      int        i3 = bMat.row_index   [index+2];
      int        j3 = bMat.column_index[index+2];
      double value3 = bMat.sp_ele      [index+2] * (*scalar);
      if (i3!=j3) {
	retMat.de_ele[i3+retMat.nCol*j3] += value3;
	retMat.de_ele[j3+retMat.nCol*i3] += value3;
      } else {
	retMat.de_ele[i3+retMat.nCol*i3] += value3;
      }
      int        i4 = bMat.row_index   [index+3];
      int        j4 = bMat.column_index[index+3];
      double value4 = bMat.sp_ele      [index+3] * (*scalar);
      if (i4!=j4) {
	retMat.de_ele[i4+retMat.nCol*j4] += value4;
	retMat.de_ele[j4+retMat.nCol*i4] += value4;
      } else {
	retMat.de_ele[i4+retMat.nCol*i4] += value4;
      }
    } // end of 'for index'
    #endif
    break;
  case rSparseMatrix::DENSE:
    if (retMat.De_Di!=rDenseMatrix::DENSE
	|| aMat.De_Di!=rDenseMatrix::DENSE) {
      rError("plus :: different matrix type");
    }
    length = retMat.nRow*retMat.nCol;
    daxpy(&length,scalar,bMat.de_ele,&IONE,retMat.de_ele,&IONE);
    break;
  case rSparseMatrix::DIAGONAL:
    if (retMat.De_Di!=rDenseMatrix::DIAGONAL
	|| aMat.De_Di!=rDenseMatrix::DIAGONAL) {
      rError("plus :: different matrix type");
    }
    daxpy(&retMat.nCol,scalar,bMat.di_ele,&IONE,
	   retMat.di_ele,&IONE);
    break;
  } // end of switch
  return _SUCCESS;
}

bool rAl::plus(rBlockVector& retVec,
	       rBlockVector& aVec,
	       rBlockVector& bVec, double* scalar)
{
  if (retVec.nBlock!=aVec.nBlock || retVec.nBlock!=bVec.nBlock) {
    rError("plus:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retVec.nBlock; ++l) {
    bool judge = plus(retVec.ele[l],aVec.ele[l],
		      bVec.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::plus(rBlockDenseMatrix& retMat,
	       rBlockDenseMatrix& aMat,
	       rBlockDenseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("plus:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = plus(retMat.ele[l],aMat.ele[l],
		      bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::plus(rBlockDenseMatrix& retMat,
	       rBlockSparseMatrix& aMat,
	       rBlockDenseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("plus:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = plus(retMat.ele[l],aMat.ele[l],
		      bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

bool rAl::plus(rBlockDenseMatrix& retMat,
	       rBlockDenseMatrix& aMat,
	       rBlockSparseMatrix& bMat,
	       double* scalar)
{
  if (retMat.nBlock!=aMat.nBlock || retMat.nBlock!=bMat.nBlock) {
    rError("plus:: different nBlock size");
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<retMat.nBlock; ++l) {
    bool judge = plus(retMat.ele[l],aMat.ele[l],
		      bMat.ele[l],scalar);
    if (judge == FAILURE) {
      total_judge = FAILURE;
    }
  }
  return total_judge;
}

// ret = a '*' (*scalar)
bool rAl::let(rVector& retVec, const char eq,
	      rVector& aVec, const char op,
	      double* scalar)
{
  switch (op) {
  case '*':
    return multiply(retVec,aVec,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '*' (*scalar)
bool rAl::let(rBlockVector& retVec, const char eq,
	      rBlockVector& aVec, const char op,
	      double* scalar)
{
  switch (op) {
  case '*':
    return multiply(retVec,aVec,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '*' (*scalar)
bool rAl::let(rBlockDenseMatrix& retMat, const char eq,
	      rBlockDenseMatrix& aMat, const char op,
	      double* scalar)
{
  switch (op) {
  case '*':
    return multiply(retMat,aMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' b*(*scalar)
bool rAl::let(rVector& retVec, const char eq,
	      rVector& aVec, const char op,
	      rVector& bVec, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retVec,aVec,bVec,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retVec,aVec,bVec,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' '*' 't' 'T' b*(*scalar)
bool rAl::let(rDenseMatrix& retMat, const char eq,
	      rDenseMatrix& aMat, const char op,
	      rDenseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  case 't':
    // ret = aMat**T * bMat
    return tran_multiply(retMat,aMat,bMat,scalar);
    break;
  case 'T':
    // ret = aMat * bMat**T 
    return multiply_tran(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool rAl::let(rDenseMatrix& retMat, const char eq,
	      rSparseMatrix& aMat, const char op,
	      rDenseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool rAl::let(rDenseMatrix& retMat, const char eq,
	      rDenseMatrix& aMat, const char op,
	      rSparseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}


// ret = a '+' '-' '*' 't' 'T' b*(*scalar)
bool rAl::let(rBlockDenseMatrix& retMat, const char eq,
	      rBlockDenseMatrix& aMat, const char op,
	      rBlockDenseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  case 't':
    // ret = aMat**T * bMat
    return tran_multiply(retMat,aMat,bMat,scalar);
    break;
  case 'T':
    // ret = aMat * bMat**T
    return multiply_tran(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool rAl::let(rBlockDenseMatrix& retMat, const char eq,
	      rBlockSparseMatrix& aMat, const char op,
	      rBlockDenseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = a '+' '-' '*' b*(*scalar)
bool rAl::let(rBlockDenseMatrix& retMat, const char eq,
	      rBlockDenseMatrix& aMat, const char op,
	      rBlockSparseMatrix& bMat, double* scalar)
{
  double minus_scalar;
  switch (op) {
  case '+':
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '-':
    if (scalar) {
      minus_scalar = -(*scalar);
      scalar = &minus_scalar;
    } else {
      scalar = &DMONE;
    }
    return plus(retMat,aMat,bMat,scalar);
    break;
  case '*':
    return multiply(retMat,aMat,bMat,scalar);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = aMat '*' '/' bVec
bool rAl::let(rVector& rVec, const char eq,
	      rDenseMatrix& aMat, const char op,
	      rVector& bVec)
{
  switch (op) {
  case '*':
    return multiply(rVec,aMat,bVec,NULL);
    break;
  case '/':
    // ret = aMat^{-1} * bVec;
    // aMat is positive definite
    // and already colesky factorized.
    return solveSystems(rVec,aMat,bVec);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rVector& aVec, const char op,
	      rVector& bVec)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aVec,bVec);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rDenseMatrix& aMat, const char op,
	      rDenseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aMat,bMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rDenseMatrix& aMat, const char op,
	      rSparseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,bMat,aMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rSparseMatrix& aMat, const char op,
	      rDenseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aMat,bMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rBlockVector& aVec, const char op,
	      rBlockVector& bVec)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aVec,bVec);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rBlockDenseMatrix& aMat, const char op,
	      rBlockDenseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aMat,bMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}
  
// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rBlockSparseMatrix& aMat, const char op,
	      rBlockDenseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,aMat,bMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

// ret = inner_product(a,b) // op = '.'
bool rAl::let(double& ret, const char eq,
	      rBlockDenseMatrix& aMat, const char op,
	      rBlockSparseMatrix& bMat)
{
  switch (op) {
  case '.':
    return getInnerProduct(ret,bMat,aMat);
    break;
  default:
    rError("let:: operator error");
    break;
  }
  return FAILURE;
}

